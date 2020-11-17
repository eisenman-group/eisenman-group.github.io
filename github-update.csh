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
# can get a DOI directly from figshare.
# can also update DOI for github site, which is
#   doi:10.5281/zenodo.3628744
# by going to https://github.com/eisenman-group/eisenman-group.github.io/releases
# can see new DOI at 
# https://zenodo.org/account/settings/github/repository/eisenman-group/eisenman-group.github.io

git add --all
git commit -m "update"
git push -u origin master
