#!/usr/bin/env python3
import os
import re

def fix_underlines(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Define the pattern for headings and underlines
    # This matches a line of text followed by a line of the same character
    pattern = r'([^\n]+)\n([=\-~\^]+)\n'
    
    def fix_match(match):
        title = match.group(1)
        underline_char = match.group(2)[0]  # Get the character used for underlining
        return f"{title}\n{underline_char * len(title)}\n"
    
    # Replace all headings with fixed underlines
    fixed_content = re.sub(pattern, fix_match, content)
    
    if content != fixed_content:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(fixed_content)
        print(f"Fixed: {file_path}")
        return True
    return False

def process_directory(directory):
    fixed_count = 0
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.rst'):
                file_path = os.path.join(root, file)
                if fix_underlines(file_path):
                    fixed_count += 1
    return fixed_count

if __name__ == "__main__":
    # Process all .rst files in the docs directory and its subdirectories
    docs_dir = os.path.dirname(os.path.abspath(__file__))
    fixed_count = process_directory(docs_dir)
    print(f"Fixed underlines in {fixed_count} files.") 