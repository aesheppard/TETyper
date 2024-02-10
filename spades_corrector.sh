#!/bin/bash
pkg_path=$(conda info | grep -m 1 'active env location' | awk '{print $NF}')
current_path=$(pwd)
file_path="$pkg_path/share/spades-3.11.1*/share/spades/pyyaml3"
cd $file_path
old_string="if not isinstance(key, collections.Hashable):"
new_string="if not isinstance(key, collections.abc.Hashable):"
sed -i "s/$old_string/$new_string/g" constructor.py
cd "$current_path"
echo "Spades has been corrrected. You can now run TETyper as normal."
echo "Spades has been corrected successfully. Do not delete this file." > correction.txt
chmod 644 correction.txt