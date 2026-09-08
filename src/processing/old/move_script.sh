#!/bin/bash

# Define source and destination directories
SOURCE_DIR="t4-data"
DEST_DIR="filtered-data"

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

echo "Copying subdirectories from '$SOURCE_DIR' to '$DEST_DIR'..."
echo "--------------------------------------------------------"

# Iterate through items in SOURCE_DIR (without trailing slashes)
for dir in "$SOURCE_DIR"/*; do
    # Ensure it is actually a directory (skips regular files inside t1-data)
    [ -d "$dir" ] || continue

    # Extract just the folder name (e.g., "15-0.4")
    folder_name=$(basename "$dir")

    # Check if the folder does NOT exist in the destination
    if [ ! -d "$DEST_DIR/$folder_name" ]; then
        echo "[COPYING] $folder_name"
        cp -r "$dir" "$DEST_DIR/"
    else
        echo "[SKIPPED] $folder_name (already exists in destination)"
    fi
done

echo "--------------------------------------------------------"
echo "Done!"