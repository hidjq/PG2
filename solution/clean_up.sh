find ./ -name '*.pickle' -type f -print -exec rm {} \;  
find ./ -name '*.json' -type f  -print -exec rm {} \;  
find ./ -name '*.tree' -type f  -print -exec rm {} \;  
rm -rf data/

echo 'Done.'