for dir in src/chapter_0*; 
    do
    CHAPNAME=$(echo $dir | cut -d '/' -f 2)
    cp $dir/figs/*.png docs/$CHAPNAME/
    done