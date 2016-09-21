
# Download different version of CASA for image testing.

cd $HOME

# CASA 4.2.2
wget https://svn.cv.nrao.edu/casa/distro/linux/old/casapy-42.2.30986-pipe-1-64b.tar.gz
tar -zxf casapy-42.2.30986-pipe-1-64b.tar.gz
rm casapy-42.2.30986-pipe-1-64b.tar.gz

casa-4.2.2="$HOME/casapy-42.2.30986-pipe-1-64b/bin/casa"

# CASA 4.3.1
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.3.1-el6.tar.gz
tar -zxf casa-release-4.3.1-el6.tar.gz
rm casa-release-4.3.1-el6.tar.gz

casa-4.3.1="$HOME/casa-release-4.3.1-el6/bin/casa"

# CASA 4.4
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.4.0-el6.tar.gz
tar -zxf casa-release-4.4.0-el6.tar.gz
rm casa-release-4.4.0-el6.tar.gz

casa-4.4="$HOME/casa-release-4.4.0-el6/bin/casa"

# CASA 4.5.3
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.5.3-el6.tar.gz
tar -zxf casa-release-4.5.3-el6.tar.gz
rm casa-release-4.5.3-el6.tar.gz

casa-4.5.3="$HOME/casa-release-4.5.3-el6/bin/casa"

# CASA 4.6
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.6.0-el6.tar.gz
tar -zxf casa-release-4.6.0-el6.tar.gz
rm casa-release-4.6.0-el6.tar.gz

casa="$HOME/casa-release-4.6.0-el6/bin/casa"
