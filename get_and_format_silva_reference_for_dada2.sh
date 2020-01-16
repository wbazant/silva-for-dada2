#!/usr/bin/bash
# Author: wbazant
# Get the SILVA release and format it as a DADA2 reference
# TODO: a special case for E. coli to match DADA2's special handling

name=$1
version=$2

if [ ! "$name" -o ! "$version" ] ; then
  echo "Usage: bash -x $0 <Silva release name/number, e.g. 138> <Silva version, e.g. SSURef_NR99>"
  exit 1
fi
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR=$SCRIPT_DIR/${name}_${version}

file_name=SILVA_${name}_${version}_tax_silva.fasta
url=https://www.arb-silva.de/fileadmin/silva_databases/release_${name}/Exports/${file_name}.gz
tax_silva=$DIR/fromProvider/${file_name}.gz

for d in fromProvider workspace final; do
  mkdir -p $DIR/$d
done

echo "$0 $1 $2" > $DIR/workspace/commands.txt

[ -s $tax_silva ] || wget --quiet $url -O ${tax_silva}

# lineage up to genus level - for assignTaxonomy in DADA2
# Skip lineages that are more than six levels (unfortunately dada2 will do strange things to them)

[ -s $DIR/final/training_set.${name}_${version}.fa.gz ] || \
  zcat $tax_silva \
  | perl -pe 's/>(.*?) (.*?);([^;]*)$/>$1\t$2\t$3/' \
  | perl -E '
$/=">"; while(<>){
  chomp;
  my ($h,@seq) = split "\n", $_; 
  next unless $h and @seq;
  my $seq = join "", @seq; 
  $seq =~ s/U/T/g;
  my ($id, $lineage, $species) = split "\t", $h;
  next if split(";", $lineage) > 6;

  $lineage =~ s/;uncultured//s;
  
  say ">$lineage;\n$seq";
}
' | gzip -c > $DIR/final/training_set.${name}_${version}.fa.gz

# species - for assignSpecies
# dada2.assignSpecies needs to match up on genus later, so we only keep species where the genus matches scientific name
# also borrow filtering from dada2.makeTaxonomyFasta_Silva and add a bit more
[ -s $DIR/final/species_assignment.${name}_${version}.fa.gz ] || \
  zcat $tax_silva \
  | perl -pe 's/>(.*?) (.*?);([^;]*)$/>$1\t$2\t$3/' \
  | perl -E '
$/=">"; while(<>){
  chomp;
  my ($h,@seq) = split "\n", $_; 
  next unless $h and @seq;
  my $seq = join "", @seq; 
  $seq =~ s/U/T/g;
  my ($id, $lineage, $species) = split "\t", $h;
  next if split(";", $lineage) > 6;
  next if $lineage =~ /uncultured$/;

  my ($genus) = reverse split (";", $lineage);
  next if $species !~ /$genus/;
  next if $species =~ /unidentified/;
  next if $species =~ /sp. /;
  next if $species =~ /sp.$/;
  next if species =~ /uncultured/;
  next if $species =~ /endosymbiont/;
  
  say ">$id $species\n$seq";
}
' | gzip -c > $DIR/final/species_assignment.${name}_${version}.fa.gz 

