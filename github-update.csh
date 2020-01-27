#!/bin/bash                                                                                                  
# update eisenman-group.github.io
# ian, 2019
#
# initally cloned repository with
#   cd ~/other/website
#   git clone https://github.com/eisenman-group/eisenman-group.github.io code
# then edited files, then ran this from ~/other/website/code/.
#
# also need to update index.html (ediff index.html ../code.html and change links).
#
# can also update doi for github site, which is
#   doi:10.5281/zenodo.3628744
# by going to https://github.com/eisenman-group/eisenman-group.github.io/releases

git add --all
git commit -m "update"
git push -u origin master
