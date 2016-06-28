#!/usr/bin/env_perl


 my ($filepath);
 
##SOFTWARE lists:RepeatScout(trf,nseg),RepeatMasker(rmblast,repeatmaskerlibraries),blast,clustalw2

# chdir "../";
# if (-e  "softwares.tar")
# {
  # system("tar zxvf softwares.tar \n");
# }
# if (-e  "softwares.tar.gz")
# {
  # system("tar zxvf softwares.tar.gz \n");
# }
# chdir "./softwares";

# my $softwaresdir=`pwd`;
# $softwaresdir=~s/\s+$|\n$|\t$//g;
# my $scriptdir=`pwd`;
# $scriptdir=~s/\s+$|\n$|\t$//g;

chdir "../";
my $tmpath=`pwd`;
$tmpath=~s/\s+$|\n$|\t$//g;

my $softwaresdir="$tmpath\/softwares";
if ( -d $softwaresdir ) {
    print "$softwaresdir exists!\n";
}
else {
    mkdir $softwaresdir;
	#$softwaresdir=$temp_dir;
    print "mkdir $softwaresdir\n";
}
###find script directory :

my $scriptdir="$tmpath\/scripts";
if ( -d $scriptdir ) {
    print "$scriptdir exists!\n";
}
else {
    mkdir $scriptdir;
    print "mkdir $scriptdir\n";
}

chdir "$softwaresdir/";
#system("sh autodownload.sh");

###write into env variable
writenvbefore();
#
####clustalw安装与配置
chdir "$softwaresdir/";
if (-e  "clustalw-2.1.tar.gz")
{;}
else
{
	system("wget http://www.clustal.org/download/current/clustalw-2.1.tar.gz &>>install.log");
}

system("tar zxvf clustalw-2.1.tar.gz \n");
chdir "$softwaresdir/clustalw-2.1";
system("./configure \n");
system("make \n");
#system(" su \n");
system("make install \n");

##安装RMBlast过程:
chdir "$softwaresdir/";
if (-e  "blast-2.2.26-x64-linux.tar.gz")
{;}
else
{system("wget http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/release/2.2.26/blast-2.2.26-x64-linux.tar.gz &>>install.log");}

if (-e  "ncbi-blast-2.2.28+-x64-linux.tar.gz")
{;}
else
{system("wget http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/blast%2b/2.2.28/ncbi-blast-2.2.28+-x64-linux.tar.gz &>>install.log");}

if (-e  "ncbi-rmblastn-2.2.28-x64-linux.tar.gz")
{;}
else
{system("wget http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/rmblast/2.2.28/ncbi-rmblastn-2.2.28-x64-linux.tar.gz &>>install.log");}

system("tar zxvf blast-2.2.26-x64-linux.tar.gz \n");
system("tar zxvf ncbi-blast-2.2.28+-x64-linux.tar.gz \n");
system("tar zxvf ncbi-rmblastn-2.2.28-x64-linux.tar.gz \n");
system("cp -R ncbi-rmblastn-2.2.28/* blast-2.2.26/ \n");
system("cp -R ncbi-blast-2.2.28+/* blast-2.2.26/ \n");
system("rm -rf ncbi-rmblastn-2.2.28 \n");
system("rm -rf ncbi-blast-2.2.28+ \n");
system("mv blast-2.2.26 ncbi-rmblastn-2.2.28 \n");
print "install RMBlast success! \n";

##RepeatMask 安装与配置
chdir "$softwaresdir/";

#system("wget http://www.repeatmasker.org/RepeatMasker-open-4-0-3.tar.gz &>>install.log");
if (-e  "RepeatMasker-open-4-0-3.tar.gz*")
{;}
else
{
	system("wget --user-agent=\"Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.2.3) Gecko/20100401 Firefox/3.6.3 (.NET CLR 3.5.30729)\" http://www.repeatmasker.org/RepeatMasker-open-4-0-3.tar.gz &>>install.log");
}
my $tmpRMf="RepeatMasker-open-4-0-3.tar.gz";
	 system("mv RepeatMasker-open-4-0-3.tar.gz* RepeatMasker-open-4-0-3.tar.gz");


system("tar xzvf RepeatMasker-open-4-0-3.tar.gz \n");
##安装 RepeatMasker Libraries
system("cp repeatmaskerlibraries-20140131.tar.gz $softwaresdir/RepeatMasker/ \n");
chdir "$softwaresdir/RepeatMasker/";
system("tar xzvf ./repeatmaskerlibraries-20140131.tar.gz \n");
##system("rm repeatmaskerlibraries-20140131.tar \n");
##运行配置脚本文件
#chdir "$softwaresdir/RepeatMasker/";
#system("perl $softwaresdir/RepeatMasker/configure \n");
print "Only left RepeatMasker to be installed in the last! \n";

##三、RepeatScout的安装与配置
#system("cd $softwaresdir \n");
chdir "$softwaresdir/";

#system("wget http://bix.ucsd.edu/repeatscout/RepeatScout-1.0.5.tar.gz &>>install.log");
if (-e  "RepeatScout-1.0.5.tar.gz")
{;}
else
{system("wget --user-agent=\"Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.2.3) Gecko/20100401 Firefox/3.6.3 (.NET CLR 3.5.30729)\" http://bix.ucsd.edu/repeatscout/RepeatScout-1.0.5.tar.gz &>>install.log");}

# if (-e "RepeatScout-1.0.5.tar.gz*")
# {
	 system("mv RepeatScout-1.0.5.tar.gz* RepeatScout-1.0.5.tar.gz");
# }

system("tar xvfz RepeatScout-1.0.5.tar.gz \n");
chdir "$softwaresdir/RepeatScout-1";
system("make \n");
print "install RepeatScout success! \n";

#write into env variable
writenvafter();

##update the CRISPRdigger script
updatescript();

	chdir "$softwaresdir/";
	my $confile="other_config.txt";
	open(CF, ">$confile")||die "can't open this $confile\n";

   chdir ;
   my $usrootdir= `pwd`;
   ##find the $RepeatMaskerdir
   # my $RepeatMaskerdir=`pwd`;
    $usrootdir=~s/\s+$|\n$|\t$//g;
    chdir "$usrootdir/bin/RepeatMasker/";
    
    my $lastip=();
    $lastip="\n ###############  The next step: Configure RepeatMasker !   ######################### \n";
    $lastip.=" You should execute the following steps(command lines) in your command windows:  \n";
    $lastip.=" 1. cd $usrootdir/bin/RepeatMasker/ \n";
    $lastip.=" 2. perl ./configure \n";
    print $lastip;
	print CF $lastip;
	
    my $confevent=();
     $confevent="\n  Congratulations! There are some remaining steps to configure the RepeatMasker !!!!!!  \n";
     $confevent.="TRF path: $usrootdir/bin/trf \n";
     $confevent.="RMBlast path: $usrootdir/bin/ncbi-rmblastn-2.2.28/bin\n";
     $confevent.="When you see the sentence 'Add a Search Engine: ......'  \n  You should select '2'  (2. RMBlast)\n";
     
     
     $confevent.="\n ############ Running an example ##############\n";
     $confevent.="You shouled modify your  directory into the script directory: \n";
     $confevent.="1. cd $scriptdir/ \n";
     $confevent.="2. perl CRISPRdigger.pl -i example.fna \n \n";   
     print $confevent;
	 print CF $confevent;
	 close(CF);
##
sub updatescript
{
    chdir ;
    my $tempdir= `pwd`;
   ##find the bash_profile
   # my $writenvdir=`pwd`;
   $tempdir=~s/\s+$|\n$|\t$//g;
    print "\n \n environment variation path:: $tempdir/bin !!! \n";
    
    chdir "$scriptdir/";
	print "scriptdir path:: $scriptdir \n";
    ##system("cd ../scripts \n");
    my $crispr="CRISPRdigger.pl";
    open(CF, "$crispr")||die "can't open this $crispr\n";
    my $tmpcrispr="tmpdigger.pl";
    open(OF, ">$tmpcrispr")||die "can't open this $tmpcrispr\n";
    
   while (<CF>)
   {
      my $eachline=$_;
      if ($eachline=~/use lib/)
      {
        s/use lib.*;/use lib "$tempdir\/perl5lib\/lib";/;
      }
      elsif ($eachline=~/filter-stage-1\.prl >/)
      {
        s/\|.*RepeatScout.*\//\| $tempdir\/bin\/RepeatScout-1\//;
      }    
      $eachline=$_;
      print OF $eachline;
    }
   close CF;
   close OF;
   
   system("rm -f $crispr \n");
   rename ("$tmpcrispr","$crispr");
    system("chmod -R 777 $crispr \n");
}

##write environment var into profile
sub writenvbefore
{
    chdir ;
   my $writenvdir= `pwd`;
   ##find the bash_profile
  # my $writenvdir=`pwd`;
  $writenvdir=~s/\s+$|\n$|\t$//g;
    ##softwares env vars
   system("cp -rf $softwaresdir/nseg/ $writenvdir/bin/ \n");
   system("cp -f $softwaresdir/trf $writenvdir/bin/ \n");
   system("cp -rf $softwaresdir/perl5lib $writenvdir/ \n");
   
   my $bashfile="bash_profile";
  
   $filepath=findfilepath("./",$bashfile);
   
   system("chmod -R 777 $filepath \n");
   open(BF, ">>$filepath")||die "can't open this $filepath\n";
   while (<BF>)
   {;}
   ##bioperl env var 
   my $biovar1="export PERL5LIB=\${PERL5LIB}:$writenvdir/perl5lib \n";
   my $biovar2="export MANPATH=$writenvdir/perl5lib \n";
   print BF $biovar1;
   print BF $biovar2; 
   ##my $softsvar="export PATH=\$PATH:$softwaresdir/ncbi-rmblastn-2.2.28/bin:$softwaresdir/RepeatMasker:$softwaresdir/nseg:$softwaresdir/RepeatScout-1:$softwaresdir/trf:$softwaresdir/ \n";
   my $softsvar="export PATH=\$PATH:$writenvdir/bin/nseg:$writenvdir/bin/trf:$writenvdir/bin/:$writenvdir/bin/ncbi-rmblastn-2.2.28/bin:$writenvdir/bin/RepeatMasker:$writenvdir/bin/RepeatScout-1 \n";
   print BF $softsvar;
   close(BF);
   system("source $filepath \n");
}

sub writenvafter
{ 
    chdir ;
   my $writenvdir= `pwd`;
   ##find the bash_profile
  # my $writenvdir=`pwd`;
  $writenvdir=~s/\s+$|\n$|\t$//g;
  ##execute authority
   system("chmod -R 777 $writenvdir/bin \n");
  
   ##softwares env vars
   system("cp -rf $softwaresdir/RepeatMasker/  $writenvdir/bin/ \n");
   system("cp -rf $softwaresdir/RepeatScout-1/  $writenvdir/bin/ \n");
   system("cp -rf $softwaresdir/ncbi-rmblastn-2.2.28/ $writenvdir/bin/ \n");
   ##modify filter-1 stage
   system("cp -f $softwaresdir/filter-stage-1.prl  $writenvdir/bin/RepeatScout-1/ \n");
   ##maybe need copy clustalw2 binary into softwares dir
   system("cp -f $softwaresdir/clustalw-2.1/src/clustalw2 $writenvdir/bin/ \n");
   
   system("rm -rf $softwaresdir/RepeatMasker/ \n");
   system("rm -rf $softwaresdir/RepeatScout-1/ \n");
   system("rm -rf $softwaresdir/ncbi-rmblastn-2.2.28/ \n");
   system("rm -rf $softwaresdir/clustalw-2.1/ \n");

   system("source $filepath \n");
}

##find the some formate file according to some part name in some dir
sub findfilepath
{
  my ($dir,$fname)=@_;
  my ($filepath);
  
  opendir(DIR,"$dir"|| die "can't open this $dir");
  local @files =readdir(DIR);
  closedir(DIR);
  for $file (@files)
  { 
    next if($file=~m/\.$/ || $file =~m/\.\.$/);
    if ($file =~/\.*?($fname)/i)
    {
       $filepath="$dir/$file";
       return $filepath;
    }   
  }
  print "Can not find the file: $fname \n";
}
