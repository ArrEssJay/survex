#!/usr/bin/perl -w
use strict;

print "*Dir <Obey\$Dir>\n";

for (@ARGV) {
   s!.*/!!;
   my $f = $_;
   if (s!^(.*)\.([chs])$!$2.$1!) {
      $f =~ s!\.!/!;
#      $f = substr($f, 0, 10);
      print "*Rename $f $_\n";
   }
}
