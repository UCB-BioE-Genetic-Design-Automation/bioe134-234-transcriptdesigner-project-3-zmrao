import csv
import random

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS
    while optimizing for codon usage and checker constraints.
    """

    def __init__(self):
        self.rbsChooser = None
        self.codon_frequencies = {}
        self.aa_to_codons = {}
        self.aa_to_weights = {}
        self.best_cai_codon = {}
        self.codon_to_aa = {}
        self.stop_codon = "TAA"

        self.max_attempts = 200

        self.forbidden_checker = None
        self.promoter_checker = None
        self.codon_checker = None
        self.internal_rbs_checker = None

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self._load_codon_usage()

        self.forbidden_checker = ForbiddenSequenceChecker()
        self.forbidden_checker.initiate()

        self.promoter_checker = PromoterChecker()
        self.promoter_checker.initiate()

        self.codon_checker = CodonChecker()
        self.codon_checker.initiate()

        self.internal_rbs_checker = InternalRBSChecker()
        self.internal_rbs_checker.initiate()

    def _load_codon_usage(self) -> None:
        codon_usage_file = "genedesign/data/codon_usage.txt"
        self.codon_frequencies = {}
        self.aa_to_codons = {}
        self.aa_to_weights = {}
        self.best_cai_codon = {}
        self.codon_to_aa = {}
        stop_codons = []

        with open(codon_usage_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 3:
                    continue

                codon = row[0].strip().upper()
                aa = row[1].strip().upper()
                usage_freq = float(row[2].strip())

                self.codon_frequencies[codon] = usage_freq
                self.codon_to_aa[codon] = aa

                if aa == "*":
                    stop_codons.append((codon, usage_freq))
                    continue

                self.aa_to_codons.setdefault(aa, []).append(codon)

        for aa, codons in self.aa_to_codons.items():
            weights = [self.codon_frequencies.get(c, 0.01) for c in codons]
            self.aa_to_weights[aa] = weights
            self.best_cai_codon[aa] = codons[weights.index(max(weights))]

        if stop_codons:
            self.stop_codon = max(stop_codons, key=lambda x: x[1])[0]

    def _choose_codon(self, aa: str) -> str:
        if aa not in self.aa_to_codons:
            raise ValueError(f"Unknown amino acid '{aa}' in peptide.")

        codons = self.aa_to_codons[aa]
        weights = self.aa_to_weights[aa]
        return random.choices(codons, weights=weights, k=1)[0]

    def _find_homopolymer_issue(self, cds: str) -> int | None:
        """
        Return a nucleotide index near a problematic homopolymer/repeat, or None.
        """
        # Long homopolymers
        for base in ("A", "T", "G", "C"):
            motif = base * 6
            idx = cds.find(motif)
            if idx != -1:
                return idx + 2

        # Repeated same codon runs like AAAAAA... from repeating the same codon
        codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
        for i in range(len(codons) - 2):
            if codons[i] == codons[i + 1] == codons[i + 2]:
                return i * 3 + 1

        return None

    def _synonymous_alternatives(self, codon: str) -> list[str]:
        aa = self.codon_to_aa.get(codon)
        if aa is None or aa == "*":
            return [codon]

        synonyms = list(self.aa_to_codons.get(aa, [codon]))
        synonyms.sort(key=lambda c: self.codon_frequencies.get(c, 0.0), reverse=True)
        return synonyms

    def _replace_codon_at_nt_index(self, codons: list[str], nt_index: int) -> list[str]:
        """
        Replace the codon covering nt_index with the best alternative synonym
        that is different from the current codon.
        """
        codon_index = nt_index // 3
        if codon_index < 0 or codon_index >= len(codons) - 1:
            return codons  # do not mutate stop codon or out of range

        old_codon = codons[codon_index]
        choices = [c for c in self._synonymous_alternatives(old_codon) if c != old_codon]
        if not choices:
            return codons

        # Try a few best alternatives first
        for new_codon in choices[:3]:
            new_codons = list(codons)
            new_codons[codon_index] = new_codon
            return new_codons

        return codons

    def _repair_candidate(self, codons: list[str]) -> list[str]:
        """
        Cheap local repair pass:
        1. Break homopolymers / repeated-codon runs
        2. Break forbidden motifs if found
        Performs only a small number of local synonymous edits.
        """
        repaired = list(codons)

        for _ in range(3):
            cds = "".join(repaired)
            changed = False

            # First: fix obvious homopolymers / repeated codons
            issue_idx = self._find_homopolymer_issue(cds)
            if issue_idx is not None:
                new_codons = self._replace_codon_at_nt_index(repaired, issue_idx)
                if new_codons != repaired:
                    repaired = new_codons
                    changed = True
                    continue

            # Second: fix forbidden motif if present
            passed_forbidden, forbidden_site = self.forbidden_checker.run(cds)
            if not passed_forbidden and forbidden_site:
                hit_idx = cds.find(forbidden_site)
                if hit_idx != -1:
                    new_codons = self._replace_codon_at_nt_index(repaired, hit_idx + len(forbidden_site) // 2)
                    if new_codons != repaired:
                        repaired = new_codons
                        changed = True
                        continue

            if not changed:
                break

        return repaired

    def _score_candidate(self, cds: str, codons: list[str], rbs) -> tuple[bool, tuple]:
        transcript_dna = rbs.utr.upper() + cds

        passed_forbidden, _ = self.forbidden_checker.run(transcript_dna)
        passed_promoter, _ = self.promoter_checker.run(transcript_dna)
        passed_internal_rbs, _ = self.internal_rbs_checker.run(cds)
        codons_above_board, codon_diversity, rare_codon_count, cai_value = self.codon_checker.run(codons)

        # Keep hairpin last because it tends to be one of the more expensive checks.
        passed_hairpin, _ = hairpin_checker(transcript_dna)

        score = (
            0 if passed_forbidden else 1,
            0 if passed_promoter else 1,
            0 if passed_hairpin else 1,
            0 if passed_internal_rbs else 1,
            0 if codons_above_board else 1,
            -cai_value,
            -codon_diversity,
            rare_codon_count,
        )

        passed_all = score[:5] == (0, 0, 0, 0, 0)
        return passed_all, score

    def run(self, peptide: str, ignores: set) -> Transcript:
        best_candidate = None
        best_score = None

        for _ in range(self.max_attempts):
            codons = [self._choose_codon(aa) for aa in peptide]
            codons.append(self.stop_codon)

            codons = self._repair_candidate(codons)
            cds = "".join(codons)

            selected_rbs = self.rbsChooser.run(cds, ignores)
            passed_all, score = self._score_candidate(cds, codons, selected_rbs)

            if passed_all:
                return Transcript(selected_rbs, peptide, codons)

            if best_score is None or score < best_score:
                best_score = score
                best_candidate = (selected_rbs, codons)

        if best_candidate is not None:
            selected_rbs, codons = best_candidate
            return Transcript(selected_rbs, peptide, codons)

        # Deterministic fallback
        codons = [self.best_cai_codon[aa] for aa in peptide]
        codons.append(self.stop_codon)
        codons = self._repair_candidate(codons)
        cds = "".join(codons)
        selected_rbs = self.rbsChooser.run(cds, ignores)
        return Transcript(selected_rbs, peptide, codons)


if __name__ == "__main__":
    peptide = "MYPFIRTARMTV"

    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)

    print(transcript)