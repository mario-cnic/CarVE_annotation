#!/usr/bin/env perl
# Output of UTRAnnotatorPGP::format_utr_effects on a synthetic multi-effect variant.
# Needs the VEP container (BaseVepPlugin); tests/perl/run_in_vep_container.sh runs it under several hash seeds.
use strict;
use warnings;
use UTRAnnotatorPGP;

my $GAIN      = '5_prime_UTR_premature_start_codon_gain_variant';
my $STOP_LOSS = '5_prime_UTR_uORF_stop_codon_loss_variant';
my $FRAME     = '5_prime_UTR_uORF_frameshift_variant';

my %flags = (uAUG_gained => $GAIN, uSTOP_lost => $STOP_LOSS, uFrameshift => $FRAME, uAUG_lost => '', uSTOP_gained => undef);
my %annotations = (
  uAUG_gained => {
    1 => {KozakContext => 'ACCATGT', KozakStrength => 'Strong', DistanceToCDS => 20, type => 'OutOfFrame_oORF', DistanceToStop => 'NA', CapDistanceToStart => 5},
    2 => {KozakContext => 'TTTATGA', KozakStrength => 'Weak', DistanceToCDS => 41, type => 'uORF', DistanceToStop => 12, CapDistanceToStart => 3},
  },
  uSTOP_lost => {1 => {AltStop => 'False', AltStopDistanceToCDS => 'NA', FrameWithCDS => 'outOfFrame', KozakContext => 'GCCATGG', KozakStrength => 'Moderate', Evidence => 'False'}},
  uFrameshift => {1 => {alt_type => 'OutOfFrame_oORF', alt_type_length => 15, ref_type => 'uORF', ref_type_length => 45, KozakContext => 'ATGATGA', KozakStrength => 'Moderate', Evidence => 'False', ref_StartDistanceToCDS => 421}},
  uAUG_lost => {}, uSTOP_gained => {},
);

my ($consequences, $entries) = UTRAnnotatorPGP->format_utr_effects(\%flags, \%annotations);
my @want_consequences = ($GAIN, $STOP_LOSS, $FRAME);
my $got = join(",", @$consequences);
die "consequences: got $got\n" unless $got eq join(",", @want_consequences);
die "entries: got " . scalar(@$entries) . ", want 4\n" unless @$entries == 4;
my $want_first = "consequence=$GAIN:CapDistanceToStart=5:DistanceToCDS=20:DistanceToStop=NA:KozakContext=ACCATGT:KozakStrength=Strong:type=OutOfFrame_oORF";
die "first entry: got $entries->[0]\n" unless $entries->[0] eq $want_first;
die "second entry is not uORF 2 of the gain\n" unless $entries->[1] =~ /^consequence=\Q$GAIN\E:.*KozakContext=TTTATGA/;
die "third entry is not the stop loss\n" unless $entries->[2] =~ /^consequence=\Q$STOP_LOSS\E:/;
die "fourth entry is not the frameshift\n" unless $entries->[3] =~ /^consequence=\Q$FRAME\E:/;

my ($none_c, $none_e) = UTRAnnotatorPGP->format_utr_effects({uAUG_gained => '', uSTOP_lost => undef}, {});
die "no effect should give empty lists\n" if @$none_c || @$none_e;

print join(",", @$consequences), "\t", join(",", @$entries), "\n";
