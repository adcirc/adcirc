#!/usr/bin/env python3
"""
ADCIRC Wiki Content Importer v4

This script fetches content from a MediaWiki site, downloads images,
and converts the content to reStructuredText format for Sphinx documentation.

This version uses Pandoc as the conversion engine from MediaWiki to RST,
with minimal pre/post-processing focused on image handling.

Usage:
    python wiki_import_v4.py --page "Page_Title" --output "output.rst" [--wiki "https://wiki.url"]

Author: ADCIRC Documentation Team
"""

import os
import re
import sys
import logging
import argparse
import urllib.parse
import tempfile
import subprocess
import shutil
from pathlib import Path

import requests
from bs4 import BeautifulSoup

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('wiki_importer')

class WikiImporter:
    """Main class for importing and converting MediaWiki content to reStructuredText."""
    
    def __init__(self, wiki_url="https://wiki.adcirc.org", output_dir="docs/user_guide", 
                 img_dir="docs/_static/images", img_url_prefix="/_static/images"):
        """Initialize the WikiImporter with the wiki URL and output directories.
        
        Args:
            wiki_url (str): URL of the MediaWiki site
            output_dir (str): Directory to save RST files
            img_dir (str): Directory to save downloaded images
            img_url_prefix (str): URL prefix for images in the documentation
        """
        self.wiki_url = wiki_url
        self.output_dir = Path(output_dir)
        self.img_dir = Path(img_dir)
        self.img_url_prefix = img_url_prefix
        
        # Try to determine the API URL
        self.api_url = self._determine_api_url()
        
        # Set up a requests session with a proper User-Agent
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'ADCIRC-DocsImporter/4.0 (https://adcirc.org; documentation@adcirc.org) python-requests/2.28.1'
        })
        
        # Check if pandoc is available
        self.pandoc_available = self._check_pandoc()
        
        # Ensure directories exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.img_dir.mkdir(parents=True, exist_ok=True)
        
        # Mapping of image URLs to local paths
        self.image_map = {}
        
        # Set to track downloaded images in the current session
        self.downloaded_images = set()

    def _determine_api_url(self):
        """Determine the MediaWiki API URL based on the wiki URL."""
        # Try standard MediaWiki API location
        api_test_url = f"{self.wiki_url}/w/api.php"
        try:
            logger.info(f"Testing API URL: {api_test_url}")
            response = requests.head(api_test_url)
            if response.status_code < 400:  # 2xx or 3xx status code
                logger.info(f"Using standard API URL: {api_test_url}")
                return api_test_url
            
            # Try without the /w/ part
            api_test_url = f"{self.wiki_url}/api.php"
            logger.info(f"Testing alternate API URL: {api_test_url}")
            response = requests.head(api_test_url)
            if response.status_code < 400:
                logger.info(f"Using alternate API URL: {api_test_url}")
                return api_test_url
            
            # Default to standard location
            logger.warning(f"Could not find API URL, defaulting to: {self.wiki_url}/w/api.php")
            return f"{self.wiki_url}/w/api.php"
            
        except Exception as e:
            logger.warning(f"Error testing API URL: {e}. Defaulting to: {self.wiki_url}/w/api.php")
            return f"{self.wiki_url}/w/api.php"
            
    def _check_pandoc(self):
        """Check if pandoc is available on the system.
        
        Returns:
            bool: True if pandoc is available, False otherwise
        """
        try:
            result = subprocess.run(
                ["pandoc", "--version"], 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                check=False
            )
            if result.returncode == 0:
                version = result.stdout.decode('utf-8').splitlines()[0]
                logger.info(f"Found pandoc: {version}")
                return True
            else:
                logger.warning("Pandoc command failed, conversion may not work properly")
                return False
        except FileNotFoundError:
            logger.warning("Pandoc not found on the system, will use fallback conversion")
            return False
            
    def _is_unique_caption(self, caption, existing_captions):
        """Check if a caption is unique compared to existing ones with normalization."""
        # Normalize the current caption
        normalized_caption = caption.lower().replace('fig ', 'figure ').replace('fig. ', 'figure ')
        
        # Check both original and normalized form
        return caption not in existing_captions and normalized_caption not in existing_captions

    def fetch_page_content(self, page_title):
        """Fetch the content of a wiki page using the MediaWiki API.
        
        Args:
            page_title (str): Title of the wiki page to fetch
            
        Returns:
            dict: Page content and metadata
        """
        logger.info(f"Fetching content for page: {page_title}")
        
        # Try with standard API endpoint
        params = {
            "action": "parse",
            "page": urllib.parse.unquote(page_title),  # Ensure proper encoding
            "format": "json",
            "prop": "text|wikitext|images|categories",
            "disabletoc": "true"
        }
        
        response = self.session.get(self.api_url, params=params)
        if response.status_code != 200:
            logger.warning(f"API request failed with status {response.status_code}, trying alternatives")
            return self._try_alternative_fetch_methods(page_title)
        
        data = response.json()
        if "error" in data:
            logger.warning(f"API reported error: {data['error']}, trying alternative methods")
            return self._try_alternative_fetch_methods(page_title)
            
        if "parse" not in data:
            raise Exception(f"Unexpected API response format for page '{page_title}'")
            
        return data["parse"]
        
    def _try_alternative_fetch_methods(self, page_title):
        """Try alternative methods to fetch page content when standard API fails.
        
        Args:
            page_title (str): Title of the wiki page to fetch
            
        Returns:
            dict: Page content and metadata
        """
        # Try the edit page to get the raw wikitext
        encoded_title = urllib.parse.quote(page_title)
        edit_url = f"{self.wiki_url}/index.php?title={encoded_title}&action=edit"
        logger.info(f"Trying to fetch content from edit page: {edit_url}")
        
        try:
            response = self.session.get(edit_url)
            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')
                textarea = soup.find('textarea', {'id': 'wpTextbox1'})
                
                if textarea:
                    # Extract the raw wikitext from the textarea
                    wikitext = textarea.string or ""
                    
                    # Create a structure similar to the parse API response
                    return {
                        "title": page_title,
                        "wikitext": {"*": wikitext},
                        "text": {"*": ""},  # Empty as we have wikitext directly
                        "images": []  # Will handle images separately
                    }
        except Exception as e:
            logger.error(f"Error fetching from edit page: {e}")
        
        # Try query API as last resort
        params = {
            "action": "query",
            "titles": page_title,
            "prop": "revisions",
            "rvprop": "content",
            "rvslots": "main",
            "format": "json"
        }
        
        logger.info(f"Trying query API for page: {page_title}")
        response = self.session.get(self.api_url, params=params)
        if response.status_code == 200:
            data = response.json()
            pages = data.get("query", {}).get("pages", {})
            
            if not pages:
                raise Exception(f"Page '{page_title}' not found")
                
            page_data = next(iter(pages.values()))
            if "missing" in page_data:
                raise Exception(f"Page '{page_title}' not found")
                
            revisions = page_data.get("revisions", [])
            if not revisions:
                raise Exception(f"No revisions found for page '{page_title}'")
                
            content = revisions[0].get("slots", {}).get("main", {}).get("*", "")
            if not content:
                content = revisions[0].get("*", "")
                
            if not content:
                raise Exception(f"No content found for page '{page_title}'")
                
            # Construct a similar structure to what parse API returns
            return {
                "title": page_title,
                "wikitext": {"*": content},
                "text": {"*": ""},  # Empty as we have wikitext directly
                "images": []  # Will handle images separately
            }
        
        raise Exception(f"Failed to fetch content for page '{page_title}'")
        
    def _extract_images_from_wikitext(self, wikitext):
        """Extract image filenames from wikitext content.
        
        Args:
            wikitext (str): MediaWiki formatted text
            
        Returns:
            list: List of image filenames
        """
        # Find image references in wiki syntax
        image_patterns = [
            r'\[\[File:([^\]\|]+)',
            r'\[\[Image:([^\]\|]+)',
            r'\{\{[^\}]*?([^\}\/]+\.(jpg|jpeg|png|gif|svg))',
            r'!\[[^\]]*\]\(([^)]+\.(jpg|jpeg|png|gif|svg))\)',
            r'<img[^>]+src=["\'](.*?)["\']',  # HTML img tags
        ]
        
        images = []
        for pattern in image_patterns:
            matches = re.findall(pattern, wikitext, re.IGNORECASE)
            for match in matches:
                if isinstance(match, tuple):  # For pattern with capture groups
                    img = match[0]
                else:
                    img = match
                img = img.strip()
                if img and img not in images:
                    if '%' in img:  # Handle URL-encoded image names
                        try:
                            img = urllib.parse.unquote(img)
                        except Exception:
                            pass
                    # Clean up image names from table markup
                    img = self._clean_image_name(img)
                    images.append(img)
        
        # Also search for MediaWiki's thumb URLs
        thumb_patterns = [
            r'(/thumb/[^"\']+)',  # Thumbnail path in MediaWiki
            r'([^"\'=]+\.(jpg|jpeg|png|gif|svg)[^"\']*)',  # Any image filename
        ]
        
        for pattern in thumb_patterns:
            matches = re.findall(pattern, wikitext, re.IGNORECASE)
            for match in matches:
                if isinstance(match, tuple):
                    img = match[0]
                else:
                    img = match
                img = img.strip()
                # Only include if it looks like an image path and isn't already in our list
                if img and any(ext in img.lower() for ext in ['.jpg', '.jpeg', '.png', '.gif', '.svg']) and img not in images:
                    # Clean up image names from table markup
                    img = self._clean_image_name(img)
                    images.append(img)
        
        logger.debug(f"Extracted images: {images}")
        return images
        
    def _clean_image_name(self, filename):
        """Clean image name by removing table markup or other unwanted characters.
        
        Args:
            filename (str): Raw image filename
            
        Returns:
            str: Cleaned image filename
        """
        # Remove leading pipe characters and underscores that might come from table markup
        cleaned = re.sub(r'^[|_\s]+', '', filename)
        
        # If we still have a File: or Image: tag, extract just the filename
        if '[[File:' in cleaned or '[[Image:' in cleaned:
            file_match = re.search(r'\[\[(?:File|Image):([^\]\|]+)', cleaned)
            if file_match:
                cleaned = file_match.group(1).strip()
                
        # Handle any other malformed filename issues
        cleaned = cleaned.replace('\\', '/').strip()
        
        # Log if we made changes
        if cleaned != filename:
            logger.debug(f"Cleaned image filename: '{filename}' -> '{cleaned}'")
            
        return cleaned

    def get_image_url(self, filename):
        """Get the full URL of an image from its filename.
        
        Args:
            filename (str): Image filename
            
        Returns:
            str: Full URL of the image
        """
        logger.info(f"Getting URL for image: {filename}")
        
        # Clean filename first to handle any table markup
        filename = self._clean_image_name(filename)
        
        # Remove 'File:' prefix if it's already there to avoid duplication
        if filename.startswith("File:"):
            api_filename = filename
            clean_filename = filename[5:]
        else:
            api_filename = f"File:{filename}"
            clean_filename = filename
            
        # Create a safe filename for local path checking
        safe_name = re.sub(r'[^\w\-\.]', '_', clean_filename)
        
        # Check if the file already exists in any subdirectory of img_dir
        for root, dirs, files in os.walk(self.img_dir):
            for file in files:
                if file.lower() == safe_name.lower():
                    # Found the file locally, construct the URL
                    rel_path = os.path.relpath(os.path.join(root, file), self.img_dir)
                    local_url = f"{self.img_url_prefix}/{rel_path}"
                    logger.info(f"Found local image file: {file} -> {local_url}")
                    return local_url
        
        # If we get here, the file doesn't exist locally, so try to get it from the wiki
        # Try MediaWiki API first
        params = {
            "action": "query",
            "titles": api_filename,
            "prop": "imageinfo",
            "iiprop": "url",
            "format": "json"
        }
        
        try:
            response = self.session.get(self.api_url, params=params)
            if response.status_code == 200:
                data = response.json()
                pages = data.get("query", {}).get("pages", {})
                for page_id, page_data in pages.items():
                    # Skip missing pages
                    if "missing" in page_data:
                        continue
                        
                    imageinfo = page_data.get("imageinfo", [])
                    if imageinfo and "url" in imageinfo[0]:
                        url = imageinfo[0].get("url")
                        logger.info(f"Found image URL via API: {url}")
                        return url
        except Exception as e:
            logger.warning(f"API error getting image URL for {filename}: {e}")
        
        # If API fails, try alternative methods
        return self._try_alternative_image_methods(filename)
        
    def _try_alternative_image_methods(self, filename):
        """Try alternative methods to get image URL when API fails.
        
        Args:
            filename (str): Image filename
            
        Returns:
            str: Full URL of the image
        """
        # Clean up filename - replace spaces with underscores (MediaWiki convention)
        clean_filename = filename.replace(' ', '_')
        encoded_filename = urllib.parse.quote(clean_filename)
        
        # Create a safe filename for local path checking
        safe_name = re.sub(r'[^\w\-\.]', '_', clean_filename)
        
        # Check if the file exists locally in any subdirectory of img_dir
        for root, dirs, files in os.walk(self.img_dir):
            for file in files:
                if file.lower() == safe_name.lower() or \
                   file.lower() == clean_filename.lower() or \
                   file.lower() == filename.lower():
                    # Found the file locally, construct the URL
                    rel_path = os.path.relpath(os.path.join(root, file), self.img_dir)
                    local_url = f"{self.img_url_prefix}/{rel_path}"
                    logger.info(f"Found local image file via alternative method: {filename} -> {local_url}")
                    return local_url
        
        # If we reach here, we didn't find the file locally, so try internet methods
        
        # Try getting image from its page
        image_page_url = f"{self.wiki_url}/File:{encoded_filename}"
        logger.info(f"Checking image page: {image_page_url}")
        
        try:
            response = self.session.get(image_page_url)
            if response.status_code == 200:
                soup = BeautifulSoup(response.text, 'html.parser')
                # Look for the full-size image - usually in a div with class 'fullImageLink'
                full_image_link = soup.select_one('.fullImageLink a') or soup.select_one('.internal a')
                if full_image_link and full_image_link.has_attr('href'):
                    url = full_image_link['href']
                    if not url.startswith(('http:', 'https:')):
                        # Handle relative URLs
                        url = urllib.parse.urljoin(self.wiki_url, url)
                    logger.info(f"Found image URL from page: {url}")
                    return url
                
                # Check for img tags with the filename
                img_tags = soup.select('img')
                for img in img_tags:
                    src = img.get('src', '')
                    if clean_filename.lower() in src.lower():
                        url = src
                        if not url.startswith(('http:', 'https:')):
                            url = urllib.parse.urljoin(self.wiki_url, url)
                        logger.info(f"Found image URL from img tag: {url}")
                        return url
        except Exception as e:
            logger.warning(f"Error checking image page: {e}")
        
        # Try several common wiki image paths
        paths_to_try = [
            # MediaWiki thumbnail paths
            f"/images/thumb/{encoded_filename[0]}/{encoded_filename[0:2]}/{encoded_filename}",
            # Direct image paths
            f"/images/{encoded_filename[0]}/{encoded_filename[0:2]}/{encoded_filename}",
            f"/images/{encoded_filename}",
            f"/w/images/{encoded_filename}",
            f"/wiki/Special:FilePath/File:{encoded_filename}"
        ]
        
        # Try with common image extensions if the file doesn't already have one
        has_extension = any(clean_filename.lower().endswith(ext) for ext in 
                           ['.png', '.jpg', '.jpeg', '.gif', '.svg', '.webp'])
        
        if not has_extension:
            extensions = ['.png', '.jpg', '.jpeg', '.gif', '.svg', '.webp']
            extended_paths = []
            for path in paths_to_try:
                for ext in extensions:
                    extended_paths.append(f"{path}{ext}")
            paths_to_try.extend(extended_paths)
        
        # Try all paths
        for path in paths_to_try:
            try:
                url = urllib.parse.urljoin(self.wiki_url, path)
                response = self.session.head(url)
                if response.status_code < 400:
                    logger.info(f"Found image URL through direct path: {url}")
                    return url
            except Exception:
                continue
                
        # If all else fails, construct a likely URL (may not work)
        possible_url = f"{self.wiki_url}/images/{encoded_filename[0]}/{encoded_filename[0:2]}/{encoded_filename}"
        logger.warning(f"Using constructed URL (may not work): {possible_url}")
        return possible_url
        
    def download_images(self, page_data, page_title, subdir=None):
        """Download images from the wiki page and return a mapping of image URLs to local paths.
        
        Args:
            page_data (dict): Page data from the API
            page_title (str): Title of the wiki page
            subdir (str, optional): Subdirectory for images
            
        Returns:
            dict: Mapping of image URLs to local paths
        """
        # Get images from page data or extract from content
        images = page_data.get("images", [])
        if not images and "text" in page_data:
            if isinstance(page_data["text"], dict) and "*" in page_data["text"]:
                images = self._extract_images_from_wikitext(page_data["text"]["*"])
            elif isinstance(page_data["text"], str):
                images = self._extract_images_from_wikitext(page_data["text"])
                
        # Also check wikitext if available
        if not images and "wikitext" in page_data:
            if isinstance(page_data["wikitext"], dict) and "*" in page_data["wikitext"]:
                images = self._extract_images_from_wikitext(page_data["wikitext"]["*"])
            elif isinstance(page_data["wikitext"], str):
                images = self._extract_images_from_wikitext(page_data["wikitext"])
        
        # Extract image URLs from HTML content as well
        if "text" in page_data and isinstance(page_data["text"], dict) and "*" in page_data["text"]:
            content = page_data["text"]["*"]
            if content.strip().startswith('<'):
                soup = BeautifulSoup(content, 'html.parser')
                for img in soup.find_all('img'):
                    src = img.get('src')
                    if src and src not in images:
                        # Clean up image names from table markup
                        src = self._clean_image_name(src)
                        images.append(src)
            
        if not images:
            logger.info("No images found in page data")
            return {}
            
        logger.info(f"Found {len(images)} images to download")
        
        # Clean up any potentially malformed image names first
        clean_images = []
        for img in images:
            clean_img = self._clean_image_name(img)
            if clean_img not in clean_images:
                clean_images.append(clean_img)
        
        if len(clean_images) != len(images):
            logger.info(f"After cleaning, found {len(clean_images)} unique images")
            images = clean_images
        
        # Create a subdirectory structure for this page within the image directory
        if subdir:
            img_subdir = self.img_dir / subdir / page_title.lower().replace(' ', '_').replace('/', '_')
            url_subdir = f"{subdir}/{page_title.lower().replace(' ', '_').replace('/', '_')}"
        else:
            img_subdir = self.img_dir / page_title.lower().replace(' ', '_').replace('/', '_')
            url_subdir = page_title.lower().replace(' ', '_').replace('/', '_')
            
        img_subdir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Images will be saved to: {img_subdir}")
        
        image_map = {}
        
        # Create placeholder paths - important for missing images
        placeholder_path = self.img_dir / "placeholder.png"
        placeholder_url = f"{self.img_url_prefix}/placeholder.png"
        
        # Create a placeholder image if it doesn't exist
        if not os.path.exists(placeholder_path):
            try:
                logger.info("Creating placeholder image")
                # Try to use PIL if available
                from PIL import Image, ImageDraw
                img = Image.new('RGB', (400, 300), color=(255, 255, 255))
                draw = ImageDraw.Draw(img)
                draw.text((150, 150), "Image Not Found", fill=(0, 0, 0))
                img.save(placeholder_path)
            except ImportError:
                # Create an empty file if PIL is not available
                with open(placeholder_path, 'w') as f:
                    f.write("")
        
        # Process all images
        for img_name in images:
            try:
                # Handle both filenames and direct URLs
                if img_name.startswith(('http:', 'https:')):
                    img_url = img_name
                    # Create a filename from the URL (handle URL encoding)
                    img_filename = urllib.parse.unquote(img_name.split('/')[-1])
                else:
                    # Clean the filename - this is crucial for MediaWiki
                    clean_name = img_name.replace(' ', '_')
                    img_filename = re.sub(r'[^\w.-]', '_', clean_name)
                    
                    # Create a safe filename
                    safe_name = re.sub(r'[^\w\-\.]', '_', img_filename)
                    if not safe_name.endswith(('.png', '.jpg', '.jpeg', '.gif', '.svg', '.webp')):
                        # Try to determine extension from content type
                        content_type = None
                        try:
                            response = self.session.head(f"{self.wiki_url}/File:{clean_name}")
                            content_type = response.headers.get('Content-Type', '').lower()
                        except Exception:
                            pass
                            
                        if 'png' in content_type:
                            safe_name += '.png'
                        elif 'jpeg' in content_type or 'jpg' in content_type:
                            safe_name += '.jpg'
                        elif 'gif' in content_type:
                            safe_name += '.gif'
                        elif 'svg' in content_type:
                            safe_name += '.svg'
                        else:
                            safe_name += '.png'  # Default to PNG
                    
                    # Determine the local path and URL
                    local_path = img_subdir / safe_name
                    local_url = f"{self.img_url_prefix}/{url_subdir}/{safe_name}"
                    
                    # Check if we've already processed this image in the current session
                    download_key = f"{img_name}:{str(local_path)}"
                    if download_key in self.downloaded_images:
                        logger.info(f"Already processed this image in current session, skipping: {img_name}")
                        # We already have the mapping, no need to add it again
                        continue
                    
                    # Check if the file already exists
                    if os.path.exists(local_path):
                        logger.info(f"Image already exists, skipping download: {local_path}")
                        # Add to the image map even though we didn't download it
                        image_map[img_name] = local_url
                        if img_name.startswith('File:'):
                            image_map[img_name[5:]] = local_url
                        
                        # Mark as processed
                        self.downloaded_images.add(download_key)
                        continue
                    
                    # If we get here, we need to get the URL and download the image
                    img_url = self.get_image_url(clean_name)
                
                if not img_url:
                    logger.warning(f"Could not determine URL for image: {img_name}")
                    # Add placeholder mapping
                    image_map[img_name] = placeholder_url
                    if img_name.startswith('File:'):
                        image_map[img_name[5:]] = placeholder_url
                    continue
                
                # Try to download the image with proper error handling
                try:
                    response = self.session.get(img_url, stream=True, timeout=10)
                    response.raise_for_status()  # Raise exception for 4xx/5xx responses
                except Exception as e:
                    logger.warning(f"Failed to download image {img_name} from {img_url}: {e}")
                    image_map[img_name] = placeholder_url
                    if img_name.startswith('File:'):
                        image_map[img_name[5:]] = placeholder_url
                    continue
                
                # Save the image to a local file
                img_ext = os.path.splitext(img_filename)[1].lower()
                if not img_ext:
                    content_type = response.headers.get('Content-Type', '').lower()
                    if 'png' in content_type:
                        img_ext = '.png'
                    elif 'jpeg' in content_type or 'jpg' in content_type:
                        img_ext = '.jpg'
                    elif 'gif' in content_type:
                        img_ext = '.gif'
                    elif 'svg' in content_type:
                        img_ext = '.svg'
                    else:
                        img_ext = '.png'  # Default to PNG
                
                # Create a safe filename
                safe_name = re.sub(r'[^\w\-\.]', '_', img_filename)
                if not safe_name.endswith(img_ext):
                    safe_name = f"{safe_name}{img_ext}"
                
                # Determine the local path and URL
                local_path = img_subdir / safe_name
                local_url = f"{self.img_url_prefix}/{url_subdir}/{safe_name}"
                
                # Save the image
                with open(local_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                
                # Add to the image map - store both by original name and by URL
                image_map[img_name] = local_url
                image_map[img_url] = local_url
                
                # Also map without File: prefix if present
                if img_name.startswith('File:'):
                    image_map[img_name[5:]] = local_url
                
                # Mark as processed
                self.downloaded_images.add(download_key)
                
                logger.info(f"Downloaded image: {img_name} -> {local_path}")
                
            except Exception as e:
                logger.error(f"Error processing image {img_name}: {e}")
                # Add placeholder for failed downloads
                image_map[img_name] = placeholder_url
                if img_name.startswith('File:'):
                    image_map[img_name[5:]] = placeholder_url
        
        return image_map
        
    def preprocess_tables(self, content):
        """Minimal preprocessor for MediaWiki tables to make them compatible with Pandoc.
        
        This specifically addresses the width attribute by moving it into the style attribute.
        Also handles colspan attributes which Pandoc doesn't support well.
        
        Args:
            content (str): MediaWiki content
            
        Returns:
            str: Preprocessed MediaWiki content
        """
        import re
        
        # Pattern to find table declarations with width attribute
        # Matches: {| attributes width = "value" more-attributes
        pattern = r'({[|](?:[^}|]*?))\s+width\s*=\s*["\']?([^"\'}\s]+)["\']?([^}|]*)'
        
        def move_width_to_style(match):
            before_width = match.group(1)  # Everything before width
            width_value = match.group(2)   # The width value
            after_width = match.group(3)   # Everything after width
            
            # Check if there's already a style attribute
            if 'style=' in before_width:
                # Add width to existing style attribute
                before_width = before_width.replace('style="', f'style="width:{width_value}; ')
                before_width = before_width.replace("style='", f"style='width:{width_value}; ")
            else:
                # Add new style attribute with width
                before_width += f' style="width:{width_value};"'
            
            return before_width + after_width
        
        # Apply the replacement
        preprocessed_content = re.sub(pattern, move_width_to_style, content)
        
        # Handle colspan in table headers
        # Convert: !.*colspan=.*| (content) to ! content !! !! !!
        colspan_pattern = r'^!.*colspan=.*\| (.*)$'
        
        def expand_colspan(match):
            cell_content = match.group(1).strip()
            return f"! {cell_content} !! !! !!"
        
        # Apply the colspan replacement - process line by line
        lines = preprocessed_content.split('\n')
        for i in range(len(lines)):
            if re.match(colspan_pattern, lines[i]):
                lines[i] = re.sub(colspan_pattern, expand_colspan, lines[i])
        
        preprocessed_content = '\n'.join(lines)
        
        logger.info("Preprocessed table width attributes and colspan for Pandoc compatibility")
        return preprocessed_content

    def preprocess_figure_tables(self, content):
        """Preprocess tables that contain figures to make them more compatible with Pandoc.
        
        This function:
        1. Identifies tables with figures
        2. Extracts all file references
        3. Extracts figure caption text (starting with "Figure" or "Fig")
        4. Replaces the table with sequential figure tags and centered caption
        
        Args:
            content (str): MediaWiki content
            
        Returns:
            str: Preprocessed MediaWiki content
        """
        import re
        
        # Define pattern to match entire tables that may contain images
        table_pattern = r'{[|](.*?)[|]}[\s\n]*'
        
        # Save original content for comparison
        original_content = content
        
        # Various patterns to find captions
        caption_patterns = [
            # Match colspan="3" align="center" or similar patterns with figure text
            r'[|]\s*colspan=["\']\d+["\']\s*align=["\']center["\']\s*[|]\s*((?:Figure|Fig)[\s\d\.]+.*?)(?:[|]|\|-|$)',
            # Match align="center" followed by figure text
            r'[|]\s*align=["\']center["\']\s*[|]\s*((?:Figure|Fig)[\s\d\.]+.*?)(?:[|]|\|-|$)',
            # Simple pattern for figure text after a cell separator
            r'[|]\s*((?:Figure|Fig)[\s\d\.]+.*?)(?:[|]|\|-|$)'
        ]
        
        tables_processed = 0
        figures_extracted = 0
        
        # Log all matches of [[File: patterns before processing
        logger.info("BEFORE PREPROCESSING - Searching for all [[File: patterns")
        file_patterns_before = re.findall(r'\[\[File:[^\]]+\]\]', content, re.IGNORECASE | re.DOTALL)
        for i, pattern in enumerate(file_patterns_before[:10]):  # Limit to first 10 for log readability
            logger.info(f"  Before pattern {i+1}: {pattern}")
        
        def process_table(match):
            nonlocal tables_processed, figures_extracted
            table_content = match.group(1)
            
            # Check if this table contains figures (File: references)
            if not '[[File:' in table_content and not '[[Image:' in table_content:
                # No figures in this table, return unchanged
                return match.group(0)
            
            tables_processed += 1
            logger.info(f"*** Processing table {tables_processed} containing File: references ***")
            logger.info(f"Table content excerpt: {table_content[:200]}...")
            
            # Log the raw table content to help debug
            logger.info(f"Raw table content: {repr(table_content[:200])}")
            
            # Extract all file references from the table - improved pattern that doesn't grab the closing brackets
            all_file_refs = []
            # New improved pattern that properly separates parts without capturing closing brackets
            file_matches = re.findall(r'\[\[(?:File|Image):([^|\]]+)(?:\|([^|\]]+))?(?:\|([^|\]]+))?(?:\|([^|\]]+))?(?:\|([^|\]]+))?(?:\|([^|\]]+))?(?:\|([^|\]]+))?\]\]', 
                                     table_content, re.IGNORECASE | re.DOTALL)
            
            logger.info(f"Found {len(file_matches)} file references in table {tables_processed}")
            
            for i, match in enumerate(file_matches):
                filename = match[0].strip()
                figures_extracted += 1
                params = [p for p in match[1:] if p]
                
                # Look for size parameter (e.g., 360px, 800px)
                size_attr = ""
                for param in params:
                    param_clean = param.strip()
                    # Remove any closing brackets that might have been captured
                    param_clean = param_clean.replace(']]', '').strip()
                    
                    if re.search(r'\d+px', param_clean):
                        size_attr = f"|{param_clean}"
                        break
                
                all_file_refs.append((filename, size_attr))
                logger.info(f"  Extracted figure {i+1}: filename='{filename}', size='{size_attr}'")
            
            # Extract any caption text using multiple patterns
            caption_text = ""
            for pattern in caption_patterns:
                caption_matches = re.findall(pattern, table_content, re.IGNORECASE | re.DOTALL)
                if caption_matches:
                    caption_text = caption_matches[0].strip()
                    logger.info(f"Found caption: {caption_text[:100]}...")
                    break
                    
            # Preserve math expressions by temporarily replacing them
            if caption_text:
                # Clean up any leading/trailing pipes or whitespace
                caption_text = caption_text.strip('| \t\n\r')
                
                # Save math expressions to restore later
                math_expressions = []
                
                def save_math(match):
                    math_expressions.append(match.group(1))
                    return f"___MATH_{len(math_expressions)-1}___"
                
                # Replace math tags with placeholders
                caption_text = re.sub(r'<math>(.*?)</math>', save_math, caption_text, flags=re.DOTALL)
                
                # Remove any remaining HTML markup
                caption_text = re.sub(r'<[^>]+>', '', caption_text)
                
                # Clean up any underscores in text (common in wiki markup for spaces)
                caption_text = re.sub(r'_+', ' ', caption_text)
                
                # Restore math expressions
                for i, expr in enumerate(math_expressions):
                    caption_text = caption_text.replace(f"___MATH_{i}___", f"<math>{expr}</math>")
                
            # Build replacement text - sequential figures followed by caption
            # Using direct MediaWiki syntax without table structure
            replacement = "\n\n"
            
            # Add each figure with proper syntax avoiding pipes inside table cells
            for i, (filename, size_attr) in enumerate(all_file_refs):
                # Use direct figure tag syntax without confusing pipe characters in table cells
                # Ensure all parts are on a single line and there are no duplicate closing brackets
                figure_tag = f"[[File:{filename.strip()}{size_attr}|center]]"
                replacement += figure_tag + "\n\n"  # Add blank line between figures to ensure proper Pandoc conversion
                logger.info(f"  Generated figure tag {i+1}: {figure_tag}")
            
            # Log that we're adding blank lines between figures to avoid them being interpreted as table cells
            logger.info("Added blank lines between figure tags to prevent Pandoc from treating them as table cells")
            
            # Add centered caption if found
            if caption_text:
                replacement += f"\n<center>{caption_text}</center>\n\n"
                logger.info(f"  Added caption: {caption_text[:100]}...")
            
            logger.info(f"  Full replacement: {repr(replacement)}")
            return replacement
        
        # Process all tables in the content
        processed_content = re.sub(table_pattern, process_table, content, flags=re.DOTALL)
        
        # Check if content was changed
        if processed_content != original_content:
            if tables_processed > 0:
                logger.info(f"Successfully preprocessed {tables_processed} tables containing {figures_extracted} figures")
                
                # Log all matches of [[File: patterns after processing
                logger.info("AFTER PREPROCESSING - Searching for all [[File: patterns")
                file_patterns_after = re.findall(r'\[\[File:[^\]]+\]\]', processed_content, re.IGNORECASE | re.DOTALL)
                for i, pattern in enumerate(file_patterns_after[:10]):  # Limit to first 10 for log readability
                    logger.info(f"  After pattern {i+1}: {pattern}")
        else:
            logger.warning("Preprocessing didn't change the content. No tables with figures were processed.")
            
        return processed_content

    def convert_with_pandoc(self, content, output_format="rst"):
        """Convert wiki content using Pandoc.
        
        Args:
            content (str): The wiki content to convert
            output_format (str): The output format (default: rst)
            
        Returns:
            str: The converted content
        """
        if not self._check_pandoc():
            logger.warning("Pandoc not available, conversion may not be optimal")
            # Return None to indicate we need to use the fallback
            return None
            
        # Apply preprocessing
        preprocessed_content = self.preprocess_tables(content)
        preprocessed_content = self.preprocess_figure_tables(preprocessed_content)
        
        # Log the MediaWiki content before passing to pandoc - focusing on image syntax
        logger.info("BEFORE PANDOC - Examining figure tags in preprocessed content")
        file_patterns = re.findall(r'\[\[File:[^\]]+\]\]', preprocessed_content, re.IGNORECASE | re.DOTALL)
        for i, pattern in enumerate(file_patterns[:10]):  # Limit to first 10 for log readability
            logger.info(f"  Figure tag {i+1} to be processed by Pandoc: {pattern}")
        
        # Log suspicious patterns that might cause problems
        suspicious_patterns = re.findall(r'[|](?:Fig|File).*?[|][^\n]*?\]\]', preprocessed_content, re.IGNORECASE | re.DOTALL)
        if suspicious_patterns:
            logger.warning(f"Found {len(suspicious_patterns)} suspicious patterns that might cause problems:")
            for i, pattern in enumerate(suspicious_patterns[:5]):
                logger.warning(f"  Suspicious pattern {i+1}: {pattern}")
                
        # Log full preprocessed content to a debug file
        try:
            with open("pandoc_input_debug.txt", "w", encoding="utf-8") as f:
                f.write(preprocessed_content)
            logger.info("Saved preprocessed content to pandoc_input_debug.txt")
        except Exception as e:
            logger.error(f"Failed to save debug file: {e}")
        
        try:
            # Create temporary files for input and output
            with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', suffix='.wiki', delete=False) as input_file:
                input_path = input_file.name
                input_file.write(preprocessed_content)
                
            output_path = f"{input_path}.{output_format}"
            
            # Build the pandoc command
            cmd = [
                "pandoc",
                "-f", "mediawiki",  # Input format
                "-t", output_format,  # Output format 
                "-o", output_path,  # Output file
                input_path          # Input file
            ]
            
            # Add additional arguments for output format
            if output_format == "rst":
                # Use list tables instead of grid tables
                cmd.extend(["--columns=80", "--standalone"])
                        
            # Run pandoc command
            logger.info(f"Running pandoc command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False  # Don't raise an exception on error
            )
            
            # Check if conversion was successful
            if result.returncode != 0:
                logger.error(f"Pandoc conversion failed: {result.stderr.decode('utf-8')}")
                return None
                
            # Read the converted content
            with open(output_path, 'r', encoding='utf-8') as f:
                converted_content = f.read()
                
            # Log the converted content to check for problems
            logger.info("AFTER PANDOC - Examining problematic patterns in converted content")
            
            # Look for malformed figure references from the examples
            malformed_refs = re.findall(r'[|]Fig.*?[|].*?\]\]', converted_content, re.IGNORECASE | re.DOTALL)
            if malformed_refs:
                logger.error(f"Found {len(malformed_refs)} malformed figure references after Pandoc conversion:")
                for i, ref in enumerate(malformed_refs[:5]):
                    logger.error(f"  Malformed ref {i+1}: {ref}")
            
            # ADD NEW LOGGING FOR TABLES - specifically looking for issues with Table 2
            logger.info("Examining table content in converted RST...")
            
            # Log all tables found in the converted content
            table_patterns = re.findall(r'(\.\.\s+table::.*?(?=\n\n\S|\Z))', converted_content, re.DOTALL)
            logger.info(f"Found {len(table_patterns)} table directives in converted content")
            
            # Look for Table 2 specifically
            table2_patterns = re.findall(r'(\.\.\s+table::.*?Table\s+2.*?(?=\n\n\S|\Z))', converted_content, re.DOTALL)
            if table2_patterns:
                logger.info("Found Table 2 content:")
                for i, table in enumerate(table2_patterns):
                    logger.info(f"Table 2 content (first 500 chars): {table[:500]}")
                    # Save Table 2 content to a separate file for detailed inspection
                    try:
                        with open("table2_debug.txt", "w", encoding="utf-8") as f:
                            f.write(table)
                        logger.info("Saved Table 2 content to table2_debug.txt")
                    except Exception as e:
                        logger.error(f"Failed to save Table 2 debug file: {e}")
            else:
                logger.warning("Table 2 not found in converted content")
            
            # Look for list-table directives that might be Table 2
            list_tables = re.findall(r'(\.\.\s+list-table::.*?(?=\n\n\S|\Z))', converted_content, re.DOTALL)
            logger.info(f"Found {len(list_tables)} list-table directives in converted content")
            
            # Check for grid tables - these might be problematic for conversion
            grid_tables = re.findall(r'(\+[-=]+\+[\s\S]*?\+[-=]+\+)', converted_content)
            logger.info(f"Found {len(grid_tables)} potential grid tables in converted content")
            if grid_tables and len(grid_tables) > 0:
                logger.info("Grid table sample (first 500 chars):")
                logger.info(grid_tables[0][:500])
            
            # Log full converted content to a debug file
            try:
                with open("pandoc_output_debug.txt", "w", encoding="utf-8") as f:
                    f.write(converted_content)
                logger.info("Saved converted content to pandoc_output_debug.txt")
            except Exception as e:
                logger.error(f"Failed to save debug file: {e}")
                
            # Clean up temporary files
            try:
                os.unlink(input_path)
                os.unlink(output_path)
            except Exception as e:
                logger.warning(f"Failed to clean up temporary files: {e}")
                
            logger.info("Pandoc conversion successful")
            return converted_content
            
        except Exception as e:
            logger.error(f"Error during pandoc conversion: {e}")
            # Return None to indicate we need to use the fallback
            return None

    def post_process_tables(self, content):
        """Convert RST grid tables within table directives to list-tables with wrap-table class attribute.
        
        This method handles the RST format where there's a blank line between the table directive
        and the grid table content, including cases where title continuation lines don't have
        proper indentation.
        
        Args:
            content (str): RST content with table directives
            
        Returns:
            str: Updated RST content with list-tables
        """
        import re
        
        # Add logging to track Table 2 specifically 
        logger.info("Beginning table post-processing, searching for tables in content...")
        
        # Look for Table 2 in the content before processing
        table2_patterns = re.findall(r'Table\s+2\..*?(?=\n\n|$)', content, re.DOTALL | re.IGNORECASE)
        if table2_patterns:
            logger.info(f"Found references to Table 2 in content before processing: {len(table2_patterns)}")
            for i, ref in enumerate(table2_patterns):
                logger.info(f"Table 2 reference {i+1}: {ref[:100]}")
                
            # Create a debug file with a section of content around Table 2
            for match in re.finditer(r'Table\s+2\..*?(?=\n\n|$)', content, re.DOTALL | re.IGNORECASE):
                start = max(0, match.start() - 500)
                end = min(len(content), match.end() + 1000)
                table2_context = content[start:end]
                try:
                    with open("table2_context_before_processing.txt", "w", encoding="utf-8") as f:
                        f.write(table2_context)
                    logger.info("Saved context around Table 2 to table2_context_before_processing.txt")
                except Exception as e:
                    logger.error(f"Failed to save Table 2 context debug file: {e}")
                break  # Just capture the first occurrence
        else:
            logger.warning("No references to Table 2 found in content before processing")
        
        # An alternative approach using a more direct string pattern matching
        # This is more robust for RST format
        tables_converted = 0
        tables_fixed = 0
        lines = content.split('\n')
        
        # Create a fresh list for the processed content
        processed_lines = []
        i = 0
        
        # Simple state machine to track table processing
        in_table = False
        table_start_idx = -1
        table_lines = []
        
        while i < len(lines):
            line = lines[i]
            
            # Check for table directive start
            if re.match(r'^\s*\.\.\s+table::', line):
                logger.info(f"Found table directive at line {i+1}: {line}")
                in_table = True
                table_start_idx = i
                table_lines = [line]  # Start collecting table lines
                
                # Check if this might be Table 2
                if "Table 2" in line or "Table 2." in line:
                    logger.info(f"This appears to be Table 2 directive: {line}")
                
                i += 1
                continue
            
            # If we're in a table, keep collecting lines
            if in_table:
                table_lines.append(line)
                
                # For Table 2, log the lines being collected
                if any("Table 2" in tl for tl in table_lines):
                    logger.debug(f"Collecting line for Table 2: {line}")
                
                # Check if we've reached the end of the table
                # Tables end at a blank line followed by text that's not part of the grid
                if not line.strip():
                    if i + 1 < len(lines):
                        next_line = lines[i + 1]
                        # If the next line is not a grid line or continuation of the table,
                        # we've reached the end of the table
                        if next_line.strip() and not (
                            re.match(r'\s*\+', next_line) or 
                            re.match(r'\s*\|', next_line) or
                            next_line.strip().startswith(':') or  # table options
                            re.match(r'\s+\S', next_line)  # indented continuation
                        ):
                            # Process the table
                            result = self._convert_table_to_list_table(table_lines)
                            
                            # Check if this is Table 2
                            is_table2 = any("Table 2" in tl for tl in table_lines)
                            
                            if result:
                                processed_lines.extend(result)
                                in_table = False
                                tables_converted += 1
                                log_msg = f"Converted table directive {tables_converted} to list-table"
                                if is_table2:
                                    log_msg += " (THIS IS TABLE 2)"
                                    logger.info(f"Table 2 converted successfully with {len(result)} lines")
                                    # Save the Table 2 conversion result for inspection
                                    try:
                                        with open("table2_conversion_result.txt", "w", encoding="utf-8") as f:
                                            f.write("\n".join(result))
                                        logger.info("Saved Table 2 conversion result to table2_conversion_result.txt")
                                    except Exception as e:
                                        logger.error(f"Failed to save Table 2 conversion result: {e}")
                                logger.info(log_msg)
                                i += 1
                                continue
                            else:
                                # If conversion fails, check if we need to fix title indentation
                                fixed_table_lines = self._fix_table_title_indentation(table_lines)
                                
                                if is_table2:
                                    logger.warning("Table 2 conversion to list-table failed")
                                    if fixed_table_lines != table_lines:
                                        logger.info("Fixed Table 2 title indentation")
                                    else:
                                        logger.warning("No title indentation fixes needed for Table 2")
                                    
                                    # Save what we have for Table 2
                                    try:
                                        with open("table2_raw_lines.txt", "w", encoding="utf-8") as f:
                                            f.write("\n".join(table_lines))
                                        logger.info("Saved raw Table 2 lines to table2_raw_lines.txt")
                                    except Exception as e:
                                        logger.error(f"Failed to save Table 2 raw lines: {e}")
                                
                                if fixed_table_lines != table_lines:
                                    tables_fixed += 1
                                    logger.info(f"Fixed table title indentation for table directive {tables_fixed}")
                                    processed_lines.extend(fixed_table_lines)
                                else:
                                    processed_lines.extend(table_lines)
                                in_table = False
                                i += 1
                                continue
                
                # Check for explicit new directive which would end the current table
                if line.strip().startswith('.. ') and not line.strip().startswith('.. table::'):
                    # Process the table before we move on
                    result = self._convert_table_to_list_table(table_lines)
                    
                    # Check if this is Table 2
                    is_table2 = any("Table 2" in tl for tl in table_lines)
                    
                    if result:
                        processed_lines.extend(result)
                        in_table = False
                        tables_converted += 1
                        log_msg = f"Converted table directive {tables_converted} to list-table"
                        if is_table2:
                            log_msg += " (THIS IS TABLE 2)"
                            logger.info(f"Table 2 converted successfully with {len(result)} lines")
                        logger.info(log_msg)
                    else:
                        # If table processing failed, check if we need to fix title indentation
                        fixed_table_lines = self._fix_table_title_indentation(table_lines)
                        
                        if is_table2:
                            logger.warning("Table 2 conversion to list-table failed (at new directive)")
                            
                        if fixed_table_lines != table_lines:
                            tables_fixed += 1
                            logger.info(f"Fixed table title indentation for table directive {tables_fixed}")
                            processed_lines.extend(fixed_table_lines)
                        else:
                            # If no changes needed, just add the original lines
                            processed_lines.extend(table_lines)
                        in_table = False
                    
                    # Add the current line (new directive) and continue
                    processed_lines.append(line)
                    i += 1
                    continue
                
                # Continue to next line while still in table
                i += 1
                continue
            
            # If not in a table, just add the line and continue
            processed_lines.append(line)
            i += 1
        
        # Handle any remaining table at the end of the file
        if in_table:
            result = self._convert_table_to_list_table(table_lines)
            
            # Check if this is Table 2
            is_table2 = any("Table 2" in tl for tl in table_lines)
            
            if result:
                processed_lines.extend(result)
                tables_converted += 1
                log_msg = f"Converted table directive {tables_converted} to list-table"
                if is_table2:
                    log_msg += " (THIS IS TABLE 2)"
                    logger.info(f"Table 2 converted successfully with {len(result)} lines")
                logger.info(log_msg)
            else:
                # If table processing failed, check if we need to fix title indentation
                fixed_table_lines = self._fix_table_title_indentation(table_lines)
                
                if is_table2:
                    logger.warning("Table 2 conversion to list-table failed (at end of file)")
                    
                if fixed_table_lines != table_lines:
                    tables_fixed += 1
                    logger.info(f"Fixed table title indentation for table directive {tables_fixed}")
                    processed_lines.extend(fixed_table_lines)
                else:
                    # If no changes needed, just add the original lines
                    processed_lines.extend(table_lines)
        
        # Look for Table 2 in the processed content
        processed_content = '\n'.join(processed_lines)
        table2_patterns_after = re.findall(r'Table\s+2\..*?(?=\n\n|$)', processed_content, re.DOTALL | re.IGNORECASE)
        if table2_patterns_after:
            logger.info(f"Found references to Table 2 in content AFTER processing: {len(table2_patterns_after)}")
            
            # Create a debug file with a section of content around Table 2 after processing
            for match in re.finditer(r'Table\s+2\..*?(?=\n\n|$)', processed_content, re.DOTALL | re.IGNORECASE):
                start = max(0, match.start() - 500)
                end = min(len(processed_content), match.end() + 1000)
                table2_context = processed_content[start:end]
                try:
                    with open("table2_context_after_processing.txt", "w", encoding="utf-8") as f:
                        f.write(table2_context)
                    logger.info("Saved context around Table 2 to table2_context_after_processing.txt")
                except Exception as e:
                    logger.error(f"Failed to save Table 2 context debug file: {e}")
                break  # Just capture the first occurrence
        else:
            logger.warning("No references to Table 2 found in content AFTER processing")
        
        if tables_converted > 0 or tables_fixed > 0:
            total_tables = tables_converted + tables_fixed
            logger.info(f"Processed {total_tables} tables: {tables_converted} converted to list-tables, {tables_fixed} had title indentation fixed")
            return processed_content
        else:
            logger.warning("No table directives were found and processed")
            return content

    def _fix_table_title_indentation(self, table_lines):
        """Fix indentation for table title continuation lines and remove bold markers.
        
        Args:
            table_lines (list): Lines of the table directive
            
        Returns:
            list: Lines with fixed indentation for title continuation
        """
        import re
        
        fixed_lines = []
        in_title = False
        directive_indent = 0
        directive_line_idx = -1
        opened_bold = False
        
        # First, identify the table directive line and capture any title content
        for i, line in enumerate(table_lines):
            if re.match(r'^\s*\.\.\s+table::', line):
                directive_line_idx = i
                in_title = True
                
                # Calculate proper indentation for continuation lines
                directive_match = re.match(r'^(\s*)\.\.\s+table::', line)
                if directive_match:
                    # Get the base indentation plus 3 spaces for directive content
                    directive_indent = len(directive_match.group(1)) + 3
                else:
                    directive_indent = 3  # Default if we can't determine
                
                # Extract and clean the title from the directive line
                title_match = re.match(r'^\s*\.\.\s+table::\s+(.*)', line)
                title = ""
                if title_match:
                    title = title_match.group(1).strip()
                
                # Check for bold markers in the title
                if title.startswith('**') and title.endswith('**'):
                    # Both opening and closing bold markers in the same line
                    title = title[2:-2]
                    logger.info(f"Removed both bold markers from table title: {title}")
                elif title.startswith('**'):
                    # Only opening bold marker - we should check for closing in continuation lines
                    title = title[2:]
                    opened_bold = True
                    logger.info(f"Removed opening bold marker from table title: {title}")
                
                # Replace the line with the modified title
                fixed_lines.append(f"{line[:line.find('::') + 2]} {title}")
                continue
            
            # If we're still processing title lines
            if in_title and i > directive_line_idx:
                # Check if this is likely a continuation of the title
                # We consider it a title continuation if it contains text and is followed by
                # either a blank line, a grid table start, or a directive option
                
                cleaned_line = line.strip()
                
                # Handle grid table start or blank line - signals end of title unless we're still looking for bold ending
                if not cleaned_line or re.match(r'\s*\+', line) or re.match(r'\s*:', line):
                    # If we previously found an opening bold marker but no closing marker,
                    # we're done with the title
                    in_title = False
                    fixed_lines.append(line)
                    continue
                
                # Check for the closing bold marker if we previously found an opening one
                if opened_bold and cleaned_line.endswith('**'):
                    cleaned_line = cleaned_line[:-2]
                    opened_bold = False
                    logger.info(f"Removed closing bold marker from continuation line: {cleaned_line}")
                
                # This is a title continuation line that needs proper indentation
                if i < len(table_lines) - 1:
                    next_line = table_lines[i + 1]
                    
                    # The line is continuation if:
                    # 1. The next line is blank, or
                    # 2. The next line starts a grid table, or
                    # 3. The next line is a directive option
                    is_continuation = (not next_line.strip() or 
                                      re.match(r'\s*\+', next_line) or 
                                      re.match(r'\s*:', next_line))
                    
                    # Definitely continuation if line not properly indented and 
                    # next line is blank, grid start, or option
                    if is_continuation and not line.startswith(' ' * directive_indent):
                        fixed_lines.append(' ' * directive_indent + cleaned_line)
                        logger.debug(f"Fixed indentation for title continuation line: '{line}' -> '{' ' * directive_indent + cleaned_line}'")
                        continue
                
                # Handle any other unindented lines that are likely title continuations
                # This is a bit more lenient to catch more cases
                if cleaned_line and not line.startswith(' ' * directive_indent):
                    # Unindented line with content that doesn't start with a special character
                    # and is right after the directive line or another continuation line
                    if not re.match(r'^\s*[.+|:]', line):
                        fixed_lines.append(' ' * directive_indent + cleaned_line)
                        logger.debug(f"Fixed indentation for additional title continuation line: '{line}' -> '{' ' * directive_indent + cleaned_line}'")
                        continue
                
                # End of title if we get here
                in_title = False
            
            # If we get here, the line is either not part of the title or already properly formatted
            fixed_lines.append(line)
        
        return fixed_lines

    def _convert_table_to_list_table(self, table_lines):
        """Helper function to convert a table directive with grid table to a list-table.
        
        Args:
            table_lines (list): Lines of the table directive and grid table
            
        Returns:
            list: Lines of the converted list-table, or None if conversion failed
        """
        import re
        
        # Check if this might be Table 2
        is_table2 = any("Table 2" in line for line in table_lines)
        if is_table2:
            logger.info("Processing Table 2 in _convert_table_to_list_table")
            # Save raw table lines for debugging
            try:
                with open("table2_raw_convert_input.txt", "w", encoding="utf-8") as f:
                    f.write("\n".join(table_lines))
                logger.info("Saved raw Table 2 convert input to table2_raw_convert_input.txt")
            except Exception as e:
                logger.error(f"Failed to save Table 2 raw convert input: {e}")
        
        # Check if we have a grid table
        has_grid_table = False
        for line in table_lines:
            if re.match(r'\s*\+[-=]+\+', line):
                has_grid_table = True
                break
        
        if not has_grid_table:
            if is_table2:
                logger.warning("Table 2 directive does not contain a grid table")
            else:
                logger.info("Table directive does not contain a grid table")
            return None
        
        # Extract the directive and title
        title = ""
        directive_line_idx = -1
        directive_indent = 0
        
        for i, line in enumerate(table_lines):
            if re.match(r'^\s*\.\.\s+table::', line):
                directive_line_idx = i
                
                # Get the title from the directive line
                title_match = re.match(r'^\s*\.\.\s+table::\s+(.*)', line)
                if title_match:
                    title = title_match.group(1).strip()
                    if is_table2:
                        logger.info(f"Table 2 title extracted: {title}")
                
                # Calculate indentation for checking continuation lines
                directive_match = re.match(r'^(\s*)\.\.\s+table::', line)
                if directive_match:
                    directive_indent = len(directive_match.group(1)) + 3
                else:
                    directive_indent = 3
                
                # Check for unindented continuation lines (look ahead until grid table start)
                j = i + 1
                while j < len(table_lines):
                    next_line = table_lines[j]
                    # Stop at grid table start or blank line followed by grid
                    if re.match(r'\s*\+', next_line) or not next_line.strip():
                        if j + 1 < len(table_lines) and re.match(r'\s*\+', table_lines[j+1]):
                            break
                    
                    # Check if this line is likely a title continuation
                    # Properly indented continuations are part of the title
                    if next_line.startswith(' ' * directive_indent) and next_line.strip():
                        title += " " + next_line.strip()
                        if is_table2:
                            logger.info(f"Added indented continuation to Table 2 title: {next_line.strip()}")
                    # Unindented lines that are likely title continuations
                    elif next_line.strip() and not next_line.startswith(' '):
                        title += " " + next_line.strip()
                        if is_table2:
                            logger.info(f"Added unindented continuation to Table 2 title: {next_line.strip()}")
                    j += 1
                
                break
        
        if is_table2:
            logger.info(f"Full Table 2 title after processing: {title}")
        
        # Find the grid table
        grid_start = -1
        for i, line in enumerate(table_lines):
            if re.match(r'\s*\+[-=]+\+', line):
                grid_start = i
                if is_table2:
                    logger.info(f"Found grid table start for Table 2 at line {i}: {line}")
                break
        
        if grid_start < 0:
            if is_table2:
                logger.warning("Could not find grid table start for Table 2")
            else:
                logger.warning("Could not find grid table start")
            return None
        
        # Extract all grid table lines
        grid_lines = []
        i = grid_start
        while i < len(table_lines):
            # Include grid lines and content lines
            if (re.match(r'\s*\+[-=]+\+', table_lines[i]) or 
                re.match(r'\s*\|', table_lines[i]) or 
                not table_lines[i].strip()):  # Include blank lines within the table
                grid_lines.append(table_lines[i])
                i += 1
            else:
                # Stop at non-grid, non-blank lines that aren't part of the table
                break
        
        if is_table2:
            logger.info(f"Extracted {len(grid_lines)} grid lines for Table 2")
            # Log the first few lines of the grid
            for i, line in enumerate(grid_lines[:5]):
                logger.info(f"Grid line {i}: {line}")
            
            # Save grid lines for debugging
            try:
                with open("table2_grid_lines.txt", "w", encoding="utf-8") as f:
                    f.write("\n".join(grid_lines))
                logger.info("Saved Table 2 grid lines to table2_grid_lines.txt")
            except Exception as e:
                logger.error(f"Failed to save Table 2 grid lines: {e}")
        
        # Now we need to process this grid table into a proper structure
        # First, identify the column positions by looking at the first separator line
        column_positions = []
        for i, line in enumerate(grid_lines):
            if re.match(r'\s*\+[-=]+\+', line):
                # Extract the positions of each + character
                for j, char in enumerate(line):
                    if char == '+':
                        column_positions.append(j)
                if is_table2 and i == 0:  # Just log the first separator line
                    logger.info(f"Table 2 column positions extracted: {column_positions}")
                break
        
        if not column_positions:
            if is_table2:
                logger.warning("Could not identify column positions in Table 2 grid table")
            else:
                logger.warning("Could not identify column positions in grid table")
            return None
        
        # Now extract all content rows from the table
        content_rows = []
        for line in grid_lines:
            if re.match(r'\s*\|', line):
                # This is a content row
                content_rows.append(line)
        
        if is_table2:
            logger.info(f"Extracted {len(content_rows)} content rows for Table 2")
            # Log a sample of content rows
            for i, row in enumerate(content_rows[:3]):
                logger.info(f"Content row {i}: {row}")
        
        # Process the content rows into a structured table
        structured_table = []
        has_header = False
        header_separator_idx = -1
        
        # First, identify if we have a header separator (line with = characters)
        for i, line in enumerate(grid_lines):
            if '=' in line and '+' in line:
                has_header = True
                header_separator_idx = i
                if is_table2:
                    logger.info(f"Found header separator for Table 2 at index {i}: {line}")
                break
        
        if is_table2:
            if has_header:
                logger.info("Table 2 has a header row")
            else:
                logger.info("Table 2 does NOT have a header row")
        
        # Process each content row into cells
        row_data = []
        current_row = []
        in_header = has_header
        header_rows = []
        normal_rows = []
        
        for i, line in enumerate(content_rows):
            # Skip empty rows
            if not line.strip():
                continue
                
            # Check if we've processed past the header separator
            if has_header and i > 0:
                header_line_idx = grid_lines.index(content_rows[0])
                if grid_lines.index(line) > header_separator_idx:
                    in_header = False
            
            # Extract cells based on column positions
            cells = []
            for j in range(len(column_positions) - 1):
                start_pos = column_positions[j] + 1
                end_pos = column_positions[j + 1]
                
                # Get the content between these positions
                if start_pos < len(line) and end_pos <= len(line):
                    cell_content = line[start_pos:end_pos].strip()
                    cells.append(cell_content)
                else:
                    cells.append("")
            
            # Log cell extraction for Table 2
            if is_table2 and i < 3:  # Just log a few rows for Table 2
                logger.info(f"Table 2 row {i} extracted cells: {cells}")
            
            # Store the cells in the appropriate section
            if in_header:
                header_rows.append(cells)
                if is_table2:
                    logger.info(f"Added cells to Table 2 header: {cells}")
            else:
                # Improved logic for handling multi-line rows
                if normal_rows:
                    # This might be a continuation of the previous row
                    
                    # Count non-empty cells in both rows
                    prev_non_empty = sum(1 for cell in normal_rows[-1] if cell.strip())
                    curr_non_empty = sum(1 for cell in cells if cell.strip())
                    
                    # Count cells where previous row is non-empty and current row is empty
                    prev_content_curr_empty = sum(1 for j in range(len(cells)) 
                                               if j < len(normal_rows[-1]) and 
                                               normal_rows[-1][j].strip() and 
                                               not cells[j].strip())
                    
                    # Count cells where current row is non-empty and previous row is empty
                    curr_content_prev_empty = sum(1 for j in range(len(cells)) 
                                               if j < len(normal_rows[-1]) and 
                                               not normal_rows[-1][j].strip() and 
                                               cells[j].strip())
                    
                    # Count cells where both rows have content
                    both_have_content = sum(1 for j in range(len(cells)) 
                                          if j < len(normal_rows[-1]) and 
                                          normal_rows[-1][j].strip() and 
                                          cells[j].strip())
                    
                    # Define criteria for a continuation row:
                    # 1. Either the current row has mostly empty cells (appears like a continuation)
                    # 2. Or there are a small number of cells where both rows have content 
                    #    (like a date spanning two lines)
                    
                    is_continuation = False
                    
                    # If most cells where the previous row has content are empty in the current row,
                    # and the current row has very few non-empty cells, it's likely a continuation
                    if prev_non_empty > 0 and curr_non_empty < len(cells) / 2:
                        is_continuation = True
                    
                    # If there's a significant number of cells where the previous row is empty
                    # but the current row has content, it's not a continuation
                    if curr_content_prev_empty > 1 and curr_content_prev_empty > both_have_content:
                        is_continuation = False
                    
                    # Special case for date spans, where both rows might have content in the same column
                    # but it's still a continuation (e.g., "08/23-08/30," in one row and "2005" in the next)
                    if both_have_content <= prev_non_empty / 2 and both_have_content <= curr_non_empty:
                        for j in range(len(cells)):
                            if j < len(normal_rows[-1]) and normal_rows[-1][j] and cells[j]:
                                # Check if it looks like a date continuation
                                if (normal_rows[-1][j].endswith(',') or 
                                    normal_rows[-1][j].endswith('-') or
                                    cells[j].strip().isdigit()):  # Check if current cell is just a year
                                    is_continuation = True
                    
                    # If the current row is mostly empty and has very little content, 
                    # it's probably a continuation
                    if curr_non_empty == 1 and prev_non_empty > 2:
                        is_continuation = True
                    
                    if is_continuation:
                        # Merge with the previous row
                        for j, cell in enumerate(cells):
                            if cell:
                                if j < len(normal_rows[-1]):
                                    if normal_rows[-1][j]:
                                        # Logic for appending content with correct spacing
                                        prev_content = normal_rows[-1][j]
                                        
                                        # Check if this looks like a continuation of a date or phrase
                                        if (prev_content.endswith(',') or 
                                            prev_content.endswith('-') or
                                            prev_content.endswith('/')):
                                            # No space needed for dates and such
                                            normal_rows[-1][j] += " " + cell
                                        else:
                                            # Add space for normal text
                                            normal_rows[-1][j] += " " + cell
                                    else:
                                        normal_rows[-1][j] = cell
                        
                        if is_table2:
                            logger.info(f"Table 2: Merged continuation row: {cells} into previous row")
                            logger.info(f"Table 2: Result after merge: {normal_rows[-1]}")
                    else:
                        normal_rows.append(cells)
                        if is_table2:
                            logger.info(f"Table 2: Added new row: {cells}")
                else:
                    normal_rows.append(cells)
                    if is_table2:
                        logger.info(f"Table 2: Added first row: {cells}")
        
        # Now process the header rows - merge them if needed
        processed_header = []
        if header_rows:
            # Initialize with the first header row
            processed_header = [cell for cell in header_rows[0]]
            
            # Merge subsequent header rows
            for row in header_rows[1:]:
                for j, cell in enumerate(row):
                    if j < len(processed_header):
                        if cell and processed_header[j]:
                            processed_header[j] += " " + cell
                        elif cell:
                            processed_header[j] = cell
            
            if is_table2:
                logger.info(f"Table 2 processed header: {processed_header}")
        
        # Log the normal rows for Table 2
        if is_table2:
            logger.info(f"Table 2 normal rows count: {len(normal_rows)}")
            for i, row in enumerate(normal_rows):
                logger.info(f"Table 2 normal row {i}: {row}")
        
        # Generate the list-table replacement
        list_table_lines = []
        
        # Clean up title - remove any beginning and ending "**" (bold markers)
        # First check for both start and end
        if title.startswith('**') and title.endswith('**'):
            title = title[2:-2]
            logger.info(f"Removed both bold markers from table title: {title}")
        # Then check for just start or end
        elif title.startswith('**'):
            title = title[2:]
            logger.info(f"Removed opening bold marker from table title: {title}")
        elif title.endswith('**'):
            title = title[:-2]
            logger.info(f"Removed closing bold marker from table title: {title}")
        
        # Start with the directive and title
        list_table_lines.append(f".. list-table:: {title}")
        list_table_lines.append("   :class: wrap-table")
        
        # Add scroll-table class for wider tables (more than 3 columns)
        num_columns = len(column_positions) - 1
        if num_columns > 3:
            # Modify the last line to add scroll-table class
            list_table_lines[-1] = list_table_lines[-1] + " scroll-table"
            if is_table2:
                logger.info(f"Added scroll-table class for wide Table 2 with {num_columns} columns")
            else:
                logger.info(f"Added scroll-table class for wide table with {num_columns} columns")
        
        # Add header-rows parameter if we have a header
        if has_header:
            list_table_lines.append("   :header-rows: 1")
        
        # Add blank line before content
        list_table_lines.append("")
        
        # Add header row if exists
        if processed_header:
            # Clean up the header cells
            clean_header = []
            for cell in processed_header:
                # Preserve bold formatting
                cell = cell.replace('**', '**')
                # Clean up any extra whitespace
                cell = ' '.join(cell.split())
                clean_header.append(cell)
            
            header_line = "   * - " + "\n     - ".join(clean_header)
            list_table_lines.append(header_line)
            
            if is_table2:
                logger.info(f"Added header line for Table 2: {header_line}")
        
        # Add data rows
        for row in normal_rows:
            # Skip rows without cell content
            if not any(cell.strip() for cell in row):
                continue
                
            # Process cells to preserve formatting
            processed_row = []
            for cell in row:
                # Preserve bold formatting
                cell = cell.replace('**', '**')
                # Fix any broken math expressions
                cell = re.sub(r'(:math:`[^`]*)\s+([^`]*`)', r'\1\2', cell)
                # Clean up any extra whitespace
                cell = ' '.join(cell.split())
                processed_row.append(cell)
            
            row_line = "   * - " + "\n     - ".join(processed_row)
            list_table_lines.append(row_line)
            
            if is_table2:
                logger.info(f"Added data row for Table 2: {processed_row}")
        
        # Add a blank line after the table
        list_table_lines.append("")
        
        if is_table2:
            logger.info(f"Final Table 2 list-table has {len(list_table_lines)} lines")
            # Save the final list table for debugging
            try:
                with open("table2_final_list_table.txt", "w", encoding="utf-8") as f:
                    f.write("\n".join(list_table_lines))
                logger.info("Saved final Table 2 list-table to table2_final_list_table.txt")
            except Exception as e:
                logger.error(f"Failed to save Table 2 final list-table: {e}")
        
        return list_table_lines

    def fix_image_references(self, content, image_map):
        """Fix image references in the converted RST content.
        
        Args:
            content (str): RST content
            image_map (dict): Mapping of image filenames to URLs
            
        Returns:
            str: Updated RST content with fixed image references
        """
        logger.info("Fixing image references in converted content")
        
        # First, check for and fix any malformed figure references that got through
        # These would be in the format |Fig_name| \|center]] - a specific pattern from our logs
        malformed_refs = re.findall(r'[|]([^|]+)[|]\s*\\[|]center\]\]', content)
        if malformed_refs:
            logger.info(f"Found {len(malformed_refs)} malformed figure references to fix")
            for filename in malformed_refs:
                # Construct a proper image directive to replace the malformed one
                img_url = None
                
                # Check if the filename is in our image map (case-insensitive)
                for img_name, url in image_map.items():
                    img_basename = os.path.basename(img_name.strip())
                    if img_basename.lower() == filename.lower() or \
                       img_basename.lower() == filename.lower().replace('_', ' '):
                        img_url = url
                        logger.info(f"Found match for malformed reference: {filename} -> {url}")
                        break
                
                if not img_url:
                    # If no match found, use placeholder
                    img_url = f"{self.img_url_prefix}/placeholder.png"
                    logger.info(f"Using placeholder for malformed reference: {filename}")
                
                # Pattern to match the malformed reference with proper regex escaping
                malformed_pattern = re.escape(f"|{filename}| \\|center]]")
                
                # Replace with a proper RST image directive
                replacement = f".. image:: {img_url}\n   :align: center\n   :alt: {filename}"
                
                # Apply the replacement
                original_content = content
                content = re.sub(malformed_pattern, replacement, content)
                
                if content != original_content:
                    logger.info(f"Fixed malformed reference: |{filename}| \\|center]] -> {replacement}")
        
        # Now continue with standard image directive processing
        # Extract all image directives
        image_directives = re.findall(r'\.\.\s+(image|figure)::\s+(.+?)$', content, re.MULTILINE)
        
        # Process each image directive
        for directive_type, path in image_directives:
            path = path.strip()
            
            # Skip if the path is already a URL or one of our paths
            if path.startswith('http') or path.startswith(self.img_url_prefix):
                continue
                
            # Try to match the path with an entry in our image map
            img_url = None
            
            # Look for exact matches by filename (strip off any directory part)
            filename = os.path.basename(path)
            
            # Create a safe filename for local path checking
            safe_name = re.sub(r'[^\w\-\.]', '_', filename)
            
            # First, check if the file exists locally in any subdirectory of img_dir
            for root, dirs, files in os.walk(self.img_dir):
                for file in files:
                    if file.lower() == safe_name.lower() or \
                       file.lower() == filename.lower() or \
                       file.lower() == filename.lower().replace('_', ' '):
                        # Found the file locally, construct the URL
                        rel_path = os.path.relpath(os.path.join(root, file), self.img_dir)
                        img_url = f"{self.img_url_prefix}/{rel_path}"
                        logger.info(f"Found local image file for reference: {filename} -> {img_url}")
                        break
                if img_url:
                    break
            
            # If not found locally, check in the image map
            if not img_url:
                # Check if the filename is in our image map (case-insensitive)
                for img_name, url in image_map.items():
                    img_basename = os.path.basename(img_name.strip())
                    if img_basename.lower() == filename.lower():
                        img_url = url
                        logger.info(f"Found image match in map: {filename} -> {url}")
                        break
            
                # Also try with File: prefix
                if not img_url:
                    for img_name, url in image_map.items():
                        if img_name.lower() == f"file:{filename}".lower() or \
                           img_name.lower() == filename.lower():
                            img_url = url
                            logger.info(f"Found image match with File: prefix: {filename} -> {url}")
                            break
            
            # If we found a match, replace the image path in the content
            if img_url:
                # Create the replacement pattern with proper regex escaping
                pattern = re.escape(f".. {directive_type}:: {path}")
                replacement = f".. {directive_type}:: {img_url}"
                
                # Apply the replacement
                content = re.sub(pattern, replacement, content)
                logger.info(f"Replaced image reference: {path} -> {img_url}")
            else:
                # If no match found, use placeholder
                placeholder_url = f"{self.img_url_prefix}/placeholder.png"
                pattern = re.escape(f".. {directive_type}:: {path}")
                replacement = f".. {directive_type}:: {placeholder_url}"
                
                # Apply the replacement
                content = re.sub(pattern, replacement, content)
                logger.info(f"Replaced missing image with placeholder: {path} -> {placeholder_url}")
        
        # Also check for direct image references outside of directives
        # This might catch some false positives but better than missing images
        for img_name, img_url in image_map.items():
            # Skip URLs in the search
            if img_name.startswith(('http:', 'https:')):
                continue
                
            # Get the filename part for matching
            img_basename = os.path.basename(img_name.strip())
            if img_basename:
                # Look for the filename in the content (outside of directives)
                # This is a bit broad but helps catch references in table cells, etc.
                content = content.replace(f" {img_basename} ", f" {img_url} ")
                
        return content

    def remove_figure_captions(self, content):
        """Remove captions from figure directives in RST content.
        
        Args:
            content (str): RST content with figure directives
            
        Returns:
            str: RST content with figure captions removed
        """
        logger.info("Removing figure captions from content")
        
        # Define pattern for figure directives with attributes and caption
        # This matches:
        # 1. The figure directive line
        # 2. Any number of attribute lines (indented, starting with :)
        # 3. A blank line
        # 4. The caption text (any non-blank lines after the blank line until the next directive or double blank line)
        figure_pattern = r'(\.\.\s+figure::\s+[^\n]+)(\n\s+:[^\n]+)*(\n\s*\n)([^\n]+(\n[^\n]+)*)'
        
        # Function to replace with just the directive and attributes
        def replace_figure(match):
            directive = match.group(1)  # The figure directive line
            attributes = match.group(2) or ''  # Any attribute lines
            
            # Return the directive and attributes, followed by a blank line
            return f"{directive}{attributes}\n"
        
        # Apply the replacement
        modified_content = re.sub(figure_pattern, replace_figure, content)
        
        # Check if any replacements were made
        if modified_content != content:
            logger.info("Successfully removed figure captions")
            return modified_content
        
        # If no replacements made with the above pattern, try a simpler pattern
        # This is a fallback for figures with different formatting
        simple_pattern = r'(\.\.\s+figure::\s+[^\n]+)(\n\s+:[^\n]+)*(\n\s*\n)([^\n]+)'
        modified_content = re.sub(simple_pattern, replace_figure, content)
        
        if modified_content != content:
            logger.info("Successfully removed figure captions with simpler pattern")
        else:
            logger.info("No figure captions found or removed")
            
        return modified_content
        
    def convert_wiki_to_rst(self, page_title, output_file=None, subdir=None):
        """Convert MediaWiki content to reStructuredText using Pandoc.
        
        Args:
            page_title (str): Title of the MediaWiki page
            output_file (str, optional): Output file path
            subdir (str, optional): Subdirectory for output files
            
        Returns:
            str: Path to the output file
        """
        logger.info(f"Converting wiki page: {page_title}")
        
        # Fetch the page content
        page_data = self.fetch_page_content(page_title)
        
        # Download images
        image_map = self.download_images(page_data, page_title, subdir)
        self.image_map = image_map
        
        # Create a placeholder image for figures that don't have corresponding images
        placeholder_path = os.path.join(self.img_dir, "placeholder.png")
        if not os.path.exists(placeholder_path):
            try:
                # Create a simple placeholder image
                logger.info("Creating placeholder image")
                try:
                    from PIL import Image, ImageDraw
                    img = Image.new('RGB', (800, 600), color=(255, 255, 255))
                    draw = ImageDraw.Draw(img)
                    draw.text((400, 300), "Image Placeholder", fill=(0, 0, 0))
                    img.save(placeholder_path)
                except ImportError:
                    # Create an empty file as fallback
                    with open(placeholder_path, 'wb') as f:
                        f.write(b'')
                logger.info(f"Created placeholder image at {placeholder_path}")
            except Exception as e:
                logger.warning(f"Failed to create placeholder image: {e}")
        
        # Extract the wikitext - prefer wikitext over HTML
        if "wikitext" in page_data and page_data["wikitext"].get("*"):
            content = page_data["wikitext"]["*"]
        elif isinstance(page_data.get("text"), dict):
            content = page_data["text"]["*"]
            # Remove HTML if we got HTML
            if content.strip().startswith('<'):
                soup = BeautifulSoup(content, 'html.parser')
                content = soup.get_text()
        else:
            content = page_data.get("text", "")
            
        # Try to convert with Pandoc first
        try:
            # Try pandoc conversion
            rst_content = self.convert_with_pandoc(content)
            
            # If pandoc conversion failed, fall back to direct conversion
            if rst_content is None:
                logger.warning("Pandoc conversion failed, falling back to direct conversion")
                
                # Try to use the v3 importer's direct conversion method
                try:
                    # Create a v3 importer instance to use its conversion method
                    # This allows us to reuse the existing direct conversion code
                    import docs._tools.wiki_import_v3 as v3
                    importer_v3 = v3.WikiImporter(
                        wiki_url=self.wiki_url,
                        output_dir=str(self.output_dir),
                        img_dir=str(self.img_dir),
                        img_url_prefix=self.img_url_prefix
                    )
                    # Use its conversion method
                    rst_content = importer_v3.process_mediawiki_to_rst(content, image_map)
                except ImportError:
                    # If unable to import v3, implement our own fallback
                    logger.warning("Could not import wiki_import_v3, using simple fallback")
                    
                    # Very simple direct conversion as fallback
                    # Remove HTML comments
                    content = re.sub(r'<!--.*?-->', '', content, flags=re.DOTALL)
                    # Convert basic headers
                    content = content.replace('== ', '\n\n').replace(' ==', '\n\n')
                    # Convert simple math
                    content = content.replace('<math>', ':math:`').replace('</math>', '`')
                    rst_content = content
            else:
                # Fix image references in Pandoc output
                rst_content = self.fix_image_references(rst_content, image_map)
                
                # Post-process tables to use list-table format with wrap-table class
                rst_content = self.post_process_tables(rst_content)
                
                # Remove figure captions
                rst_content = self.remove_figure_captions(rst_content)
                
        except Exception as e:
            logger.error(f"Error converting content to RST: {e}")
            # Create a simplified version with error message
            error_msg = f"\n\nERROR DURING CONVERSION: {str(e)}\n\n"
            rst_content = f"Raw content could not be fully processed due to conversion error:\n\n{error_msg}\n\nPartial content follows:\n\n"
            
            # Try to clean up the content a bit
            try:
                # Remove HTML comments
                content = re.sub(r'<!--.*?-->', '', content, flags=re.DOTALL)
                # Convert basic headers
                content = content.replace('== ', '\n\n').replace(' ==', '\n\n')
                # Convert simple math
                content = content.replace('<math>', ':math:`').replace('</math>', '`')
                rst_content += content
            except Exception:
                # If that fails, return the raw content
                rst_content += content
        
        # Add title and metadata
        title = page_title.replace('_', ' ')
        header = f"{title}\n{'=' * len(title)}\n\n"
        metadata = f".. meta::\n   :description: {title} in ADCIRC\n   :keywords: adcirc, {title.lower()}\n\n"
        
        # CSS for wrap-table
        wrap_table_css = """
.. raw:: html

   <style>
   .wrap-table th, .wrap-table td {
     white-space: normal !important;
     word-wrap: break-word !important;
     max-width: 100% !important;
     overflow-wrap: break-word !important;
     hyphens: auto !important;
   }
   </style>
"""
        
        # Final content with metadata, title, and CSS
        final_content = metadata + header + rst_content + "\n\n" + wrap_table_css
        
        # Determine the output file
        if not output_file:
            slug = page_title.lower().replace(' ', '_').replace('/', '_')
            output_path = os.path.join(self.output_dir, subdir if subdir else "", f"{slug}.rst")
            output_file = os.path.normpath(output_path)
        
        # Create directory if it doesn't exist
        output_dir = os.path.dirname(output_file)
        if output_dir:  # Only create directory if there's a dirname
            os.makedirs(output_dir, exist_ok=True)
        
        # Write the output file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(final_content)
            
        logger.info(f"Conversion complete. Output file: {output_file}")
        return output_file


def main():
    """Main function to parse arguments and run the importer."""
    parser = argparse.ArgumentParser(description="Import and convert MediaWiki content to reStructuredText using Pandoc.")
    
    # Required arguments
    parser.add_argument("--page", required=True, help="Title of the wiki page to import")
    
    # Optional arguments
    parser.add_argument("--output", help="Output file path")
    parser.add_argument("--wiki", default="https://wiki.adcirc.org", help="MediaWiki URL (without /w/api.php or similar)")
    parser.add_argument("--img-dir", default="_static/images", help="Image directory path")
    parser.add_argument("--img-url-prefix", default="/_static/images", help="Image URL prefix in documentation")
    parser.add_argument("--subdir", help="Subdirectory for output files and images")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    # Set log level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Set default subdir if not provided
    if not args.subdir:
        args.subdir = os.path.dirname(args.output) if args.output else None
    
    # Create a WikiImporter instance
    importer = WikiImporter(
        wiki_url=args.wiki,
        output_dir=os.path.dirname(args.output) if args.output else "docs/user_guide",
        img_dir=args.img_dir,
        img_url_prefix=args.img_url_prefix
    )
    
    # Run the conversion
    output_file = importer.convert_wiki_to_rst(
        page_title=args.page,
        subdir=args.subdir,
        output_file=args.output
    )
    
    print(f"Conversion complete. Output file: {output_file}")

if __name__ == "__main__":
    main() 