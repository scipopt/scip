#!/usr/bin/perl
#
# Script that converts a graph in GML format to DIMACS format. All
# visual data is ignored.

my $narg = @ARGV;

if ( $narg != 1 )
{
   printf("usage: <.> <GML file>\n");
   exit(1);
}


# try to open file
open FILE, $ARGV[0] or die $!;

my $nnodes = 0;
my $nedges = 0;

my $foundnode = 0;
my $foundedge = 0;
my $inedge = 0;
my $innode = 0;

# count number of edges and nodes
while(<FILE>)
{
   chomp;

   if ( $_ =~ /node/ )
   {
      $foundnode = 1;
   }

   if ( $_ =~ /edge/ )
   {
      $foundedge = 1;
   }

   if ( $_ =~ /\[/ )
   {
      if ( $foundnode == 1 )
      {
	 $foundnode = 0;
	 $innode = 1;
	 ++$nnodes;
      }

      if ( $foundedge == 1 )
      {
	 $foundedge = 0;
	 $inedge = 1;
	 ++$nedges;
      }
   }

   if ( $_ =~ /\]/ )
   {
      if ( $innode == 1 )
      {
	 $innode = 0;
      }

      if ( $inedge == 1 )
      {
	 $inedge = 0;
      }
   }
}

# close and reopen file
close(FILE);
open FILE, $ARGV[0] or die $!;

$foundnode = 0;
$foundedge = 0;
$inedge = 0;
$innode = 0;

my $nodeid = -1;
my $nodeweight = 0;
my $edgesource = -1;
my $edgetarget = -1;


# output header
printf("p edge %d %d\n", $nnodes, $nedges);

while(<FILE>)
{
   chomp;

   if ( $_ =~ /node/ )
   {
      $foundnode = 1;
   }

   if ( $_ =~ /edge/ )
   {
      $foundedge = 1;
   }

   if ( $_ =~ /\[/ )
   {
      if ( $foundnode == 1 )
      {
	 $foundnode = 0;
	 $innode = 1;
      }

      if ( $foundedge == 1 )
      {
	 $foundedge = 0;
	 $inedge = 1;
      }
   }

   # parse node
   if ( $innode == 1 && $_ =~ /id/ )
   {
      @array = split(/\s+/);
      $nodeid = $array[2];
   }

   if ( $innode == 1 && $_ =~ /weight/ )
   {
      @array = split(/\s+/);
      $nodeweight = $array[2];
   }

   # parse edge
   if ( $inedge == 1 && $_ =~ /source/ )
   {
      @array = split(/\s+/);
      $edgesource = $array[2];
   }

   if ( $inedge == 1 && $_ =~ /target / )
   {
      @array = split(/\s+/);
      $edgetarget = $array[2];
   }

   if ( $_ =~ /\]/ )
   {
      if ( $innode == 1 )
      {
	 printf("n %d %g\n", $nodeid, $nodeweight);
	 $innode = 0;

	 # for debugging:
	 $nodeid = -1;
	 $nodeweight = -1;
      }

      if ( $inedge == 1 )
      {
	 printf("e %d %d\n", $edgesource, $edgetarget);
	 $inedge = 0;

	 # for debugging:
	 $edgesource = -1;
	 $edgetarget = -1;
      }
   }
}
