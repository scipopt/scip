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
my $nedgeskip = 0;
my $maxnode = -1;

my $foundnode = 0;
my $foundedge = 0;
my $inedge = 0;
my $innode = 0;

my $nodeid = -1;
my $nodeweight = 0;
my $edgesource = -1;
my $edgetarget = -1;

# store edges and nodes in order to avoid parallel edges
my %E;
my %V;

# count number of edges and nodes
while(<FILE>)
{
   chomp;

   # found node keyword
   if ( $_ =~ /node/ )
   {
      $foundnode = 1;
   }

   # found edge keyword
   if ( $_ =~ /edge/ )
   {
      $foundedge = 1;
   }

   # found start of block
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

   # found end of block
   if ( $_ =~ /\]/ )
   {
      if ( $innode == 1 )
      {
	 $V{$nodeid} = $nodeweight;
	 if ( $nodeid > $maxnode )
	 {
	    $maxnode = $nodeid;
	 }

	 $innode = 0;
	 ++$nnodes;

	 # for debugging:
	 $nodeid = -1;
	 $nodeweight = -1;
      }

      if ( $inedge == 1 )
      {
	  # possibly swap nodes
	  if ( $edgesource > $edgetarget )
	  {
	     my $h = $edgesource;
	     $edgesource = $edgetarget;
	     $edgetarget = $h;
	  }

	 if ( ! exists $E{$edgesource}{$edgetarget} )
	 {
	    $E{$edgesource}{$edgetarget} = 1;
	    ++$nedges;
	 }
	 else
	 {
	    ++$nedgeskip;
	 }
	 $inedge = 0;

	 # for debugging:
	 $edgesource = -1;
	 $edgetarget = -1;
      }
   }

   # parse node
   if ( $innode == 1 )
   {
      # parse node id
      if ( $_ =~ /id/ )
      {
	 @array = split(/\s+/);
	 $nodeid = $array[2];
      }

      # parse node weight
      if ( $_ =~ /weight/ )
      {
	 @array = split(/\s+/);
	 $nodeweight = $array[2];
      }
   }

   # parse edge
   if ( $inedge == 1 )
   {
      # parse source node
      if ( $_ =~ /source/ )
      {
	 @array = split(/\s+/);
	 $edgesource = $array[2];
      }

      # parse target node
      if ( $_ =~ /target / )
      {
	 @array = split(/\s+/);
	 $edgetarget = $array[2];
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

# output header
printf("p edge %d %d\n", $maxnode, $nedges);
printf("c %d skipped parallel edges.\n", $nedgeskip);

# output nodes
while (my ($key, $value) = each(%V) )
{
   printf("n %d %g\n", $key, $value);
}

# output edges
while (my ($source, $value) = each(%E) )
{
   while (my ($target, $dummy) = each($E{$source}) )
   {
      printf("e %d %d\n", $source, $target);
   }
}
