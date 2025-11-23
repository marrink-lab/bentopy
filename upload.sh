#!/bin/bash
set -ue

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RESET='\033[0m'

version=$1
repository=${2:-testpypi} # Default is testpipy.

wheels="dist/*$version*manylinux*.whl"
sources="dist/*$version.tar.gz"

echo -e "Version: ${GREEN}$1${RESET}"
echo -e "Found ${BLUE}wheels${RESET}:"
ls $wheels
echo -e "Found ${BLUE}sources${RESET}:"
ls $sources

echo -e "Ready to upload to ${RED}$repository${RESET}"
read -p "Type 'yes' to continue:  " confirmation
if [[ "$confirmation" != "yes" ]]; then
    echo "Upload cancelled."
    exit 1
fi

set -x
twine upload --verbose $wheels $sources --repository=$repository
