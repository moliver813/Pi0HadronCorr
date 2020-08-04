# Detect OS
# From: http://stackoverflow.com/a/12099167
ifeq ($(shell uname -s),Darwin)
  extraLinkOption=
else
  extraLinkOption=-Wl,--no-as-needed
endif

default:  compare

compare: compare.cc
	g++ $(USER_DEFINED) -std=c++11 $(extraLinkOption) compare.cc -o $@ -I`root-config --incdir` -I$(MYINSTALL)/include `root-config --libs` 

