for dir in src/chapter_*; 
    do
    CHAPNAME=$(echo $dir | cut -d '/' -f 2)
    cp $dir/figs/*.png docs/$CHAPNAME/
    done