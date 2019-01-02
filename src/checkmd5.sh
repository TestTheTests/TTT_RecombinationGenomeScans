#!/bin/bash

for file in "$@"; do
  generatedhash=$(md5 -q "$file")
  storedhash=$(< "$file".md5)
  if [[ $generatedhash != $storedhash ]]; then
    echo "Hash for file '$file 'does not match"
  fi
done