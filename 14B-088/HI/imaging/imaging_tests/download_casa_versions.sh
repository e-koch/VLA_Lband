
# Download different version of CASA for image testing.

cd $HOME

# CASA 4.2.2
wget https://svn.cv.nrao.edu/casa/distro/linux/old/casapy-42.2.30986-pipe-1-64b.tar.gz
tar -zxf casapy-42.2.30986-pipe-1-64b.tar.gz
rm casapy-42.2.30986-pipe-1-64b.tar.gz

mv casapy-42.2.30986-pipe-1-64b/bin/casa casapy-42.2.30986-pipe-1-64b/bin/casa-4.2.2

echo 'export PATH="$PATH:/home/ekoch/casapy-42.2.30986-pipe-1-64b/bin"' >> $HOME/.bashrc

# CASA 4.3.1
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.3.1-el6.tar.gz
tar -zxf casa-release-4.3.1-el6.tar.gz
rm casa-release-4.3.1-el6.tar.gz

mv casa-release-4.3.1-el6/bin/casa casa-release-4.3.1-el6/bin/casa-4.3.1

echo 'export PATH="$PATH:/home/ekoch/casa-release-4.3.1-el6/bin"' >> $HOME/.bashrc

# CASA 4.4
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.4.0-el6.tar.gz
tar -zxf casa-release-4.4.0-el6.tar.gz
rm casa-release-4.4.0-el6.tar.gz

mv casa-release-4.4.0-el6/bin/casa casa-release-4.4.0-el6/bin/casa-4.4

echo 'export PATH="$PATH:/home/ekoch/casa-release-4.4.0-el6/bin"' >> $HOME/.bashrc

# CASA 4.5.3
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.5.3-el6.tar.gz
tar -zxf casa-release-4.5.3-el6.tar.gz
rm casa-release-4.5.3-el6.tar.gz

mv casa-release-4.5.3-el6/bin/casa casa-release-4.5.3-el6/bin/casa-4.5.3

echo 'export PATH="$PATH:/home/ekoch/casa-release-4.5.3-el6/bin"' >> $HOME/.bashrc

# CASA 4.6
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.6.0-el6.tar.gz
tar -zxf casa-release-4.6.0-el6.tar.gz
rm casa-release-4.6.0-el6.tar.gz

mv casa-release-4.6.0-el6/bin/casa casa-release-4.6.0-el6/bin/casa-4.6

echo 'export PATH="$PATH:/home/ekoch/casa-release-4.6.0-el6/bin"' >> $HOME/.bashrc

# CASA 4.7
wget https://svn.cv.nrao.edu/casa/distro/linux/release/el6/casa-release-4.7.0-el6.tar.gz
tar -zxf casa-release-4.7.0-el6.tar.gz
rm casa-release-4.7.0-el6.tar.gz

mv casa-release-4.7.0-el6/bin/casa casa-release-4.7.0-el6/bin/casa-4.7

echo 'export PATH="$PATH:/home/ekoch/casa-release-4.7.0-el6/bin"' >> $HOME/.bashrc
