#!/usr/bin/env python3

"""
Code to source the ADCIRC version from the git commit hash and tag
and optionally create a version file for use within the code
"""

from typing import Tuple


class AdcircVersion:
    def __init__(self, directory: str):
        """
        Initialize the ADCIRC version object

        Args:
            directory (str): The directory to create the version file in or the directory to search for the version file
        """
        self.__directory = directory
        version_dict = AdcircVersion.__query_adcirc_version(self.__directory)
        self.__tag = version_dict["tag"]
        self.__hash = version_dict["hash"]

    def tag(self) -> str:
        """
        Get the ADCIRC version tag

        Returns:
            str: The ADCIRC version tag
        """
        return self.__tag

    def hash(self) -> str:
        """
        Get the ADCIRC version hash

        Returns:
            str: The ADCIRC version hash
        """
        return self.__hash

    @staticmethod
    def __query_adcirc_version(directory: str) -> dict:
        """
        Code to source the ADCIRC version from the git commit hash and tag

        Args:
            directory (str): The directory to create the version file in or the directory to search for the version file

        Returns:
            None
        """
        import os

        if AdcircVersion.__in_git_repo():
            git_hash, git_tag = AdcircVersion.__get_git_hash()
        elif os.path.exists(os.path.join(directory, "version.F")):
            filename = os.path.join(directory, "version.F")
            git_hash, git_tag = AdcircVersion.__parse_version_file(filename)
        else:
            git_hash, git_tag = "unknown", "unknown"

        return {"tag": git_tag, "hash": git_hash}

    def generate_version_file(self, directory: str) -> None:
        """
        Creates the version file as version.F

        Args:
            directory (str): The directory to create the version file in

        Returns:
            None
        """
        import os

        file_name = os.path.join(directory, "version.F")
        with open(file_name, "w") as f:
            f.write("      module version\n")
            f.write(
                '      character(80), parameter:: ADC_VERSION = "{:s}"\n'.format(
                    self.__tag
                )
            )
            f.write(
                '      character(80), parameter:: ADC_HASH = "{:s}"\n'.format(
                    self.__hash
                )
            )
            f.write("      integer, save:: FileFmtMajor = 1\n")
            f.write("      integer, save:: FileFmtMinor = 2\n")
            f.write("      integer, save:: FileFmtRev = 0\n")
            f.write("      end module\n")

    @staticmethod
    def __in_git_repo() -> bool:
        """
        Check if we are in a git repository or an exported tarball

        Returns:
            bool: True if we are in a git repository, False otherwise
        """
        import subprocess

        cmd = ["git", "status"]
        try:
            subprocess.check_output(cmd, stderr=subprocess.DEVNULL)
            return True
        except subprocess.CalledProcessError:
            return False

    @staticmethod
    def __get_git_hash() -> Tuple[str, str]:
        """
        Get the git commit hash

        Returns:
            Tuple[str, str]: The commit hash and tag
        """
        import subprocess

        cmd = ["git", "describe", "--always", "--tags"]
        git_tag = subprocess.check_output(cmd).strip().decode("utf-8")

        cmd = ["git", "rev-parse", "HEAD"]
        git_hash = subprocess.check_output(cmd).strip().decode("utf-8")

        return git_hash, git_tag

    @staticmethod
    def __parse_version_file(filename: str) -> Tuple[str, str]:
        """
        Parse the existing version file for the hash and tag

        Args:
            filename (str): The version file name
        """
        with open(filename, "r") as f:
            for line in f:
                if "ADC_VERSION" in line:
                    tag = line.split("=")[1].strip().strip('"')
                elif "ADC_HASH" in line:
                    hash = line.split("=")[1].strip().strip('"')
        return hash, tag


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create the ADCIRC version file")
    parser.add_argument(
        "--create-version-file",
        help="Create the version file",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--tag", help="Print the version to standard out", action="store_true"
    )
    parser.add_argument(
        "--hash", help="Print the tag to standard out", action="store_true"
    )
    parser.add_argument(
        "--directory",
        help="The directory to create the version file in",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    if not args.directory and args.create_version_file:
        raise ValueError("Must provide a directory to create the version file in")

    versioning = AdcircVersion(args.directory)

    if args.create_version_file:
        versioning.generate_version_file(args.directory)

    if args.tag:
        print(versioning.tag())
    elif args.hash:
        print(versioning.hash())
