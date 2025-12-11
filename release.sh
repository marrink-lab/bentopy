#!/bin/bash
set -ueo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RESET='\033[0m'

convert_version() {
    echo "$1" | sed -E 's/-alpha\./a/; s/-beta\./b/; s/-rc\./rc/'
}

get_version() {
    grep '^version =' $1 | cut -d'"' -f2
}

report_bump() {
    echo -e "Bumping version in $1: ${RED}$2${RESET} -> ${GREEN}$3${RESET}"
}

check_git_tree() {
    if git diff-index --quiet HEAD; then
        echo -e "Tree is ${GREEN}clean${RESET}."
    else
        echo -e "Tree is ${RED}dirty${RESET}."
        echo    "Commit the changes first."
        exit 1
    fi
}

change_version() {
    sed -i -e "s/^version = \".*\"/version = \"$2\"/" $1
}

version=$1
py_version=$(convert_version $version)

old_version=$(get_version Cargo.toml)
old_py_version=$(get_version pyproject.toml)

# Report on the changes we're making.
report_bump Cargo.toml     $old_version    $version
report_bump pyproject.toml $old_py_version $py_version

# Check if the source tree is clean.
check_git_tree

# Set the new version.
change_version Cargo.toml     $version
change_version pyproject.toml $py_version

# Check the project and update Cargo.lock.
cargo check

# Make a version bump commit.
git add Cargo.toml Cargo.lock pyproject.toml
git commit -m "Bump version to $version"

# Remove old build files from the dist directory.
rm -f dist/*

# Build the tar ball.
python -m build
# Create compiled wheels to distribute.
sudo ./compile-release.sh

# Upload to testpypi.
./upload.sh $py_version testpypi
