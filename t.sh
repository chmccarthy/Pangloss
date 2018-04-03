#!/bin/bash

# First thisngs first, try to find what version of ruby we have
INSTALL_RUBY=0
RUBY_WEBSITE="http://ftp.ruby-lang.org/pub/ruby"
RUBY_V=`ruby -v`
# capture the exit value
EV=$?

if test $EV -eq 127; then
	INSTALL_RUBY=1
	RUBY_V="ruby-1.9.3-p0"
	RUBY_MAJMIN="1.9"
else
	# check the version here, if it's too old, we'll need to get a newer
	# version
	RUBY_MAJMIN=`echo $RUBY_V | sed "s/ruby \([0-9]\.[0-9]\).*/\1/g"`
	RUBY_V=`echo $RUBY_V | sed "s/ruby \([0-9]\.[0-9]\.[0-9]\)p\([0-9]*\).*/ruby-\1-p\2/g"`
fi

if test $INSTALL_RUBY -eq 1; then
	# download and install ruby
	echo "Installing ruby v1.9.3p0..."
	echo "Downloading $RUBY_V.tar.gz"
	curl "$RUBY_WEBSITE/$RUBY_MAJMIN/$RUBY_V.tar.gz" > $RUBY_V.tar.gz
	tar -xzf $RUBY_V.tar.gz
	cd  $RUBY_V
	# I'm assuming that you have access to the necessary development
	# tools to compile ruby from source
	./configure >/dev/null
	YAML_INST=`make 2>&1 | grep "libyaml is missing"`
	if test !-z $YAML_INST; then
		echo "LibYAML not found, installing..."
		# We need to install libyaml first
		echo "Downloading yaml-0.14.tar.gz"
		curl http://pyyaml.org/download/libyaml/yaml-0.1.4.tar.gz > yaml-0.14.tar.gz
		tar -xzf yaml-0.14.tar.gz
		cd yaml-0.1.4
		./configure >/dev/null
		make >/dev/null
		make install >/dev/null
		cd ..
		make clean >/dev/null 2>&1
		make >/dev/null 2>&1
	fi
	make install >/dev/null
	cd ..
fi

# Now, check to make sure Imagemagick is installed...
INSTALL_IM=0
IM_WEBSITE="http://www.imagemagick.org/download/ImageMagick.tar.gz"
IM_V=`convert -version`
EV=$?

if test $EV -eq 127; then
	INSTALL_IM=1
fi

if test $INSTALL_IM -eq 1; then
	echo "Installing ImageMagick..."
	echo "Downloading ImageMagick.tar.gz"
	curl $IM_WEBSITE > ImageMagick.tgz
	IM_DIR=`tar -xvzf ImageMagick.tgz 2>&1 | sed 's/\/.*.//g' | tail -1 | sed 's/.*[ \t]//g'`
	cd $IM_DIR
	./configure >/dev/null
	make >/dev/null 2>&1
	make install >/dev/null
	cd ..
fi
	

INSTALL_GEMS=0
INSTALL_RMAGICK=0
GEM_LIST=`gem list`
EV=$?

if test $EV -eq 127; then
	INSTALL_GEMS=1
fi

if test -z "`echo $GEM_LIST | grep rmagick`"; then
	INSTALL_RMAGICK=1
fi

if test $INSTALL_GEMS -eq 1; then
	echo "Installing RubyGems v1.8.12..."
	echo "Downloading rubygems-1.8.12.tgz"
	# download and install rubygems for the appropriate version here
	curl http://production.cf.rubygems.org/rubygems/rubygems-1.8.12.tgz > rubygems-1.8.12.tgz
	tar -xzf rubygems-1.8.12.tgz
	cd rubygems-1.8.12
	ruby setup.rb >/dev/null
	cd ..	
fi

if test $INSTALL_RMAGICK -eq 1; then
	echo "Installing RMagick..."
	# Use rubygems to install rmagick
	gem install rmagick >/dev/null 2>&1
fi

echo "All programs installed"
#OK, you should be good to go now!

