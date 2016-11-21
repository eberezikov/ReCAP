#!/usr/bin/perl -w

use strict;

`/home/fedot/src/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i assembled_ctgs.00.fa -o assembled_ctgs.00.nr99.fa -r 0 -c 0.99 -T 0 -M 0 2>log.01 1>&2`;

