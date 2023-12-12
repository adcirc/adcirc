#!/usr/bin/env python3

"""
Code to export the ADCIRC git repository into a tarball
"""


def export_adcirc() -> None:
    """
    Export the ADCIRC git repository into a tarball

    Returns:
        None
    """
    import os
    import tarfile
    import shutil
    import tempfile
    import adcirc_version
    import subprocess

    # get the version
    version = adcirc_version.AdcircVersion("..")

    # create a temporary directory
    temp_dir = tempfile.mkdtemp()

    # use the git archive command to export the repository and exclude the .circle directory
    os.chdir("..")
    cmd = [
        "git",
        "archive",
        "{:s}".format(version.hash()),
        "|",
        "tar",
        "-x",
        "-C",
        "{:s}".format(temp_dir),
    ]
    out = subprocess.run(" ".join(cmd), shell=True, check=True)
    os.chdir("scripts")

    # check if there was an error during export
    if out.returncode != 0:
        raise RuntimeError("Error exporting ADCIRC repository")

    # remove the .circleci directory
    shutil.rmtree(os.path.join(temp_dir, os.path.join(temp_dir, ".circleci")))

    # create the version file
    version.generate_version_file(temp_dir)

    # create the tarball
    tarball_name = "adcirc-{}.tar.gz".format(version.tag())
    with tarfile.open(tarball_name, "w:gz") as tarball:
        tarball.add(temp_dir, arcname="adcirc")

    # Generate checksums
    generate_checksums(tarball_name)

    # Print the checksums to the screen for the user
    print("ADCIRC Package Checksums:")
    for algorithm, checksum in generate_checksums(tarball_name).items():
        print("   {:s}: {:s}".format(algorithm, checksum))


def generate_checksums(tarball_name: str) -> dict:
    """
    Generate checksums for the tarball using md5, sha1, and sha256

    Args:
        tarball_name: name of the tarball

    Returns:
        Dictionary of checksums
    """
    import hashlib

    checksums = {}
    for algorithm in ["md5", "sha1", "sha256"]:
        hash_obj = hashlib.new(algorithm)
        with open(tarball_name, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_obj.update(chunk)
        checksums[algorithm] = hash_obj.hexdigest()

    return checksums


if __name__ == "__main__":
    export_adcirc()
