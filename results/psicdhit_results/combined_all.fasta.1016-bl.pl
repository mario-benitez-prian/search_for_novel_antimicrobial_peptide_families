#!/usr/bin/perl
$host = shift;
$instance = shift;
$arg = shift;

#### random sleep, rand() can be a fraction of second
select(undef,undef,undef,rand());

if ($arg) {
  @ids = split(/,/, $arg);
}
else {
  while(1) {
    if (opendir(DDIR, "../results/psicdhit_results/combined_all.fasta-seq")) { 
      @ids = grep {/^\d+$/} readdir(DDIR);
      last;
    }
    else {
      sleep(1);
    }
  }
}

foreach $id (@ids) {

  next unless (-e "../results/psicdhit_results/combined_all.fasta-seq/$id");
  next if (-e "../results/psicdhit_results/combined_all.fasta-seq/$id.lock");
  $cmd = `touch ../results/psicdhit_results/combined_all.fasta-seq/$id.lock`;

  if (50) {
    $cmd = `blastp -outfmt 6 -db ./../results/psicdhit_results/combined_all.fasta.1016 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query ../results/psicdhit_results/combined_all.fasta-seq/$id -out ../results/psicdhit_results/combined_all.fasta-bl/$id`;
    $cmd =                         `/home/mario/software/cd-hit-v4.8.1-2019-0228/psi-cd-hit/psi-cd-hit.pl -J parse_blout_multi ../results/psicdhit_results/combined_all.fasta-bl/$id -c 0.3 -ce -1 -aS 0.7 -aL 0 -G 1 -prog blastp -bs 0 >> ../results/psicdhit_results/combined_all.fasta-blm/$host.$instance`;
  }
  elsif (1) {
    $cmd = `blastp -outfmt 6 -db ./../results/psicdhit_results/combined_all.fasta.1016 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query ../results/psicdhit_results/combined_all.fasta-seq/$id | /home/mario/software/cd-hit-v4.8.1-2019-0228/psi-cd-hit/psi-cd-hit.pl -J parse_blout ../results/psicdhit_results/combined_all.fasta-bl/$id -c 0.3 -ce -1 -aS 0.7 -aL 0 -G 1 -prog blastp -bs 1`;
  }
  else {
    $cmd = `blastp -outfmt 6 -db ./../results/psicdhit_results/combined_all.fasta.1016 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query ../results/psicdhit_results/combined_all.fasta-seq/$id -out ../results/psicdhit_results/combined_all.fasta-bl/$id`;
    $cmd =                         `/home/mario/software/cd-hit-v4.8.1-2019-0228/psi-cd-hit/psi-cd-hit.pl -J parse_blout ../results/psicdhit_results/combined_all.fasta-bl/$id -c 0.3 -ce -1 -aS 0.7 -aL 0 -G 1 -prog blastp -bs 0`;
  }
  $cmd = `rm -f  ../results/psicdhit_results/combined_all.fasta-seq/$id`;
  $cmd = `rm -f  ../results/psicdhit_results/combined_all.fasta-seq/$id.lock`;
}

($tu, $ts, $cu, $cs) = times();
$tt = $tu + $ts + $cu + $cs;
$cmd = `echo $tt >> ../results/psicdhit_results/combined_all.fasta-seq/host.$host.$instance.cpu`;

