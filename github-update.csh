#!/bin/bash                                                                                                  
# update eisenman-group.github.io
# ian, 2019
#
# initally cloned repository with
#   cd ~/other/website
#   git clone https://github.com/eisenman-group/eisenman-group.github.io code
# then edited files, then ran this from ~/other/website/code/.
#
# when updating, also need to update index.html (ediff index.html ../code.html and change links).
#
# can get a DOI directly from figshare.
# can also update DOI for github site, which is
#   doi:10.5281/zenodo.3628744
# by going to https://github.com/eisenman-group/eisenman-group.github.io/releases
# can see new DOI at 
# https://zenodo.org/account/settings/github/repository/eisenman-group/eisenman-group.github.io
#
# I am using a (classic) personal access token, which is a long a password that's good for just one year. The token is
# stored in passwords.txt. Instructions for getting a new one at
# https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token#creating-a-token
#
# When things go awry, can use
#   git reset --hard origin/master

git add --all
git commit -m "update"
git push -u origin master
