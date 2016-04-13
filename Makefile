#
#    Copyright (C) 2016 University of Southern California
#                       Andrew D. Smith and Ting Chen
#
#    Authors: Haifeng Chen and Ting Chen
#
#    This file is part of the S3.
#
#    S3 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    S3 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with S3.  If not, see <http://www.gnu.org/licenses/>.
#

S3 = $(shell pwd)

all:
	@make -C src S3=$(S3) OPT=1

install:
	@make -C src S3=$(S3) OPT=1 install

test:
	@make -C src S3=$(S3) test
.PHONY: test

clean:
	@make -C src S3=$(S3) clean
.PHONY: clean

distclean: clean
	@rm -rf $(S3)/bin
	@rm -rf $(S3)/include
.PHONY: distclean
