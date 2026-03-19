class InternalRBSChecker:
    """
    Checks for potential internal ribosome binding sites (Shine-Dalgarno motifs)
    followed by a start codon within a CDS.
    """

    def __init__(self):
        self.sd_motifs = []
        self.start_codons = []
        self.min_spacer = 0
        self.max_spacer = 0

    def initiate(self) -> None:
        """
        Initializes common Shine-Dalgarno-like motifs and spacer rules.
        """
        # Common SD-like motifs for E. coli (conservative set).
        self.sd_motifs = [
            "AGGAGG",
            "GGAGG",
            "AGGA",
            "GAGG",
            "AAGGA",
            "GGAGA",
        ]
        self.start_codons = ["ATG", "GTG", "TTG"]
        self.min_spacer = 4
        self.max_spacer = 9

    def run(self, cds: str) -> tuple[bool, str | None]:
        """
        Scans the CDS for internal SD motifs with a valid spacer and start codon.

        Returns:
            (bool, str | None): True if no internal RBS found, else False and the hit.
        """
        seq = cds.upper()

        for motif in self.sd_motifs:
            mlen = len(motif)
            start = 0
            while True:
                idx = seq.find(motif, start)
                if idx == -1:
                    break
                for spacer in range(self.min_spacer, self.max_spacer + 1):
                    codon_start = idx + mlen + spacer
                    if codon_start + 3 <= len(seq):
                        codon = seq[codon_start:codon_start + 3]
                        if codon in self.start_codons:
                            hit = seq[idx:codon_start + 3]
                            return False, hit
                start = idx + 1

        return True, None


if __name__ == "__main__":
    checker = InternalRBSChecker()
    checker.initiate()
    example = "AAAGGAGGAAAAATGAAACCC"
    print(checker.run(example))
