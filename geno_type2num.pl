#!/usr/bin/perl -w
use strict;

sub allele_stat {
  my @list = @_;
  my $p = pop(@list);
  my @allelelist = @list;
  my $n = 0;
  foreach(@allelelist){
    if ($_ eq $p) {
      #print "$_";
      $n++;
    }
  }
  return $n;
}

#print "GGMU!\n";
my $data_dir = $ARGV[0];
my $snp_info = $ARGV[1];
my $sample_info = $ARGV[2];
#my $final_output = $ARGV[3];

#print "$data_dir,$snp_info,$sample_info,$final_output\n";
my @snp_ABI=();
#my %sam_stat;
my $line;
my @line_list;
my $flag = 0;
my $i=0;
my $j=0;
my $vic = 0;
my $both = 1;
my $fam = 2;
open(SNPF, "<", $snp_info) || die "can not open the snp file!\n";
open(OUT, ">", "geno_num.temp") || die "can not open the temp file!\n";

while(<SNPF>) {
  chomp;
  $line = $_;
  if ($line =~ /rs\d+/) {
    @line_list = split('\t', $line);
    push(@snp_ABI,($line_list[5], $line_list[6]));
  }
}
print OUT "Sample ID\tStat\t";

for($i=0;$i<=$#snp_ABI;$i=$i+2){
  print OUT "$snp_ABI[$i]\t";
}
print OUT "\n";

close(SNPF);
open(SAMF, "<", $sample_info) || die "can not open the sample file!\n";

while(<SAMF>){
  chomp;
  $line = $_;
  if ($line =~ /^\d+/) {
    @line_list = split('\t', $line);
    #$sam_stat{$line_list[0].'_'.$line_list[3]} = $line_list[4];
    if ($line_list[3] eq 'NTC') {
      `grep '^$line_list[2]\tN' $data_dir/AZ-*plate$line_list[1]* > temp.txt`;
    }
    else {
      `grep '^$line_list[2]\t$line_list[3]' $data_dir/AZ-*plate$line_list[1]* > temp.txt`
    }
    if ($line_list[4] eq 'SZ'){
      print OUT "$line_list[3]\t1\t";
    }
    else {
      print OUT "$line_list[3]\t0\t";
    }
    for($i=0;$i<=$#snp_ABI;$i=$i+2){
      $flag=0;
      open(TEMP, "<", "temp.txt") || die "can not open the temp file!\n";
      while(<TEMP>) {
        chomp;
        $line = $_;
        @line_list = split('\t', $line);
        if ($snp_ABI[$i+1] eq 'NA'){
          $snp_ABI[$i+1]=$snp_ABI[$i];
        }
        else {
          $snp_ABI[$i+1]=~s/^C_+//g;
        }
        if ($line_list[0] =~ /$snp_ABI[$i+1]/){
          print OUT "$line_list[5]\t";
          $flag = 1;
          last;
#          if ($line_list[5]=~/vic/) {
#            $flag=1;
#            print "0\t";
#            last;
#          }
#          elsif ($line_list[5]=~/Both/) {
#            print "1\t";
#            $flag=1;
#            last;
#          }
#          elsif ($line_list[5]=~/fam/) {
#            print "2\t";
#            $flag=1;
#            last;
#          }
#          else {
#            print "-1\t";
#            $flag=1;
#            last;
#          }
        }
      }
      if ($flag==0) {
        print OUT "-2\t";
      }
      close(TEMP);
    }
    print OUT "\n";
  }
}

close(SAMF);
close(OUT);

open(TEMP, "<", "geno_num.temp") || die "can not open the temp file!\n";

while(<TEMP>){
  chomp;
  $line = $_;
  if (/^Sample/){
    print "$line\n";
  }
  else {
    @line_list = split('\t', $line);
    print "$line_list[0]\t$line_list[1]\t";
    if (allele_stat(@line_list[2..$#line_list], 'vic')>allele_stat(@line_list[2..$#line_list], 'fam')) {
      $vic = 0;
      $both = 1;
      $fam = 2;
    }
    else {
      $vic = 2;
      $fam = 0;
      $both = 1;
    }
    foreach(@line_list[2..$#line_list]){
      if ($_ eq 'vic'){
        print "$vic\t";
      }
      elsif ($_ eq 'fam'){
        print "$fam\t";
      }
      elsif ($_ eq 'Both'){
        print "$both\t";
      }
      elsif ($_ eq '-2'){
        print "-2\t";
      }
      else {
        print "-1\t";
      }
    }
    print "\n";
    #last;
  }
}

close(TEMP);
