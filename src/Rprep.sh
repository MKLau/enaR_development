#!/usr/bin/env bash
echo 'Removing development files from'
pwd
echo 'Are you sure (y/n)?'
read prompt

if [ "$prompt" == "y" ]; then 
    find . -name '.Rhistory' -print0 | xargs -0 rm -rf
    find . -name '*~' -print0 | xargs -0 rm -rf
    find . -name '*.aux' -print0 | xargs -0 rm -rf
    find . -name '*.bbl' -print0 | xargs -0 rm -rf
    find . -name '*.blg' -print0 | xargs -0 rm -rf
    find . -name '*.log' -print0 | xargs -0 rm -rf
    find . -name '*.out' -print0 | xargs -0 rm -rf
    find . -name '*.toc' -print0 | xargs -0 rm -rf
    echo 'Done!'
    exit 1
else 
    echo 'Ok, quitting.'
    exit 1
fi
