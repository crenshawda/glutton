(ns glutton.lexicon)

; Ref: http://en.wikipedia.org/wiki/Genetic_code 'Inverse Table' portion
; Note: I have replaced uracil with thymine for ease of comparison (dna -> rna)
; Note: N is not considered here, yet, so codons containing it are moot at the moment.
; This is a hash in the form of {:<codon> :<amino-acid>, ...}
(def *codon-translation-matrix* 
     {:GCT :A :GCC :A :GCA :A :GCG :A
      :TTA :L :TTG :L :CTT :L :CTC :L :CTA :L :CTG :L
      :CGT :R :CGC :R :CGA :R :CGG :R :AGA :R :AGG :R
      :AAA :K :AAG :K
      :AAT :N :AAC :N
      :ATG :M ; This is synonymous with the :START codon
      :GAT :D :GAC :D
      :TTT :F :TTC :F
      :TGT :C :TGC :C
      :CCT :P :CCC :P :CCA :P :CCG :P
      :CAA :Q :CAG :Q
      :TCT :S :TCC :S :TCA :S :TCG :S :AGT :S :AGC :S
      :GAA :E :GAG :E
      :ACT :T :ACC :T :ACA :T :ACG :T
      :GGT :G :GGC :G :GGA :G :GGG :G
      :TGG :W
      :CAT :H :CAC :H
      :TAT :Y :TAC :Y
      :ATT :I :ATC :I :ATA :I
      :GTT :V :GTC :V :GTA :V :GTG :V
      :TAA :. :TGA :. :TAG :.} ; This is synonymous with the :STOP codon
     )

; These are characters for convenience of direct comparison
(def *nucleotide-base-pair-dictionary*
     {
      ; 3 intermolecular H bonds
      \A \T
      \T \A
      ; 2 intermolecular H bonds
      \G \C
      \C \G
      }
     )

(def *amino-acid-dictionary*
     {:. {:average-mass 0, :monoisotopic-mass 0}
      :A {:average-mass 71.0788, :monoisotopic-mass 71.03711}
      :C {:average-mass 103.1448, :monoisotopic-mass 161.01466}
      :D {:average-mass 115.0886, :monoisotopic-mass 115.02694}
      :E {:average-mass 129.1155, :monoisotopic-mass 129.04259}
      :F {:average-mass 147.1766, :monoisotopic-mass 147.06841}
      :G {:average-mass 57.052, :monoisotopic-mass 57.02146}
      :H {:average-mass 137.1412, :monoisotopic-mass 137.05891}
      :I {:average-mass 113.1595, :monoisotopic-mass 113.08406}
      :K {:average-mass 128.1742, :monoisotopic-mass 128.09496}
      :L {:average-mass 113.1595, :monoisotopic-mass 113.08406}
      :M {:average-mass 131.1986, :monoisotopic-mass 131.04049}
      :N {:average-mass 114.1039, :monoisotopic-mass 114.04293}
      :P {:average-mass 97.1167, :monoisotopic-mass 97.05276}
      :Q {:average-mass 128.1308, :monoisotopic-mass 128.05858}
      :R {:average-mass 156.1876, :monoisotopic-mass 156.10111}
      :S {:average-mass 87.0782, :monoisotopic-mass 87.03203}
      :T {:average-mass 101.1051, :monoisotopic-mass 101.04768}
      :V {:average-mass 99.1326, :monoisotopic-mass 99.06841}
      :W {:average-mass 186.2133, :monoisotopic-mass 186.07931}
      :Y {:average-mass 163.176, :monoisotopic-mass 163.06333}}
     )

(def *assorted-constants*
     {:average-amino-mass 125.44726659484284}
     )

; See: http://en.wikipedia.org/wiki/CHON
(def *chon-constants*
     {:oxygen {:monoisitopic-mass 15.99491463
	       :average-mass 15.99940494}
      :nitrogen {:monoisotopic-mass 14.00307400
		 :average-mass 14.00674309}
      :hydrogen {:monoisotopic-mass 1.00782504
		 :average-mass 1.00794076}
      :carbon {:monoisotopic-mass 12.00000000
	       :average-mass 12.01073590}
      }
     )

; This is only defined separately from the above because water is a molecule--
; Plus I couldn't immediately figure out how to define by composing previous
; values in the same hash.
; NOTE: Brian Risk- TODO this .98 is from my observations of theoretical
; and measured peptide mass differences.  check validity.
(def *water-constants*
     {:monoisotopic-mass (+
			  ((*chon-constants* :oxygen) :monoisitopic-mass)
			  (* 2 ((*chon-constants* :hydrogen) :monoisotopic-mass))) ; H2O, get it?
      :average-mass (+
		     ((*chon-constants* :oxygen):average-mass)
		     (* 2 ((*chon-constants* :hydrogen) :average-mass)))}
     )

; Def need to ask Brian how he came up with this...
; NOTE: Brian Risk- indicies in our amino acid list which define
; trypsin cleavages.
(def *trypsin-cleavages*
     {
      :no-cleavages-before [20 23 22 21]
      :cleavages [8 10 24 25 26 0 2]
      }
     )

(def *hmm-state-constants*
     {
      :number-of-ions 12
      :b-ion 0
      :y-ion 1
      :a-ion 2
      :j-ion 3
      :b17-ion 4
      :y17-ion 5
      :b18-ion 6
      :y18-ion 7
      :imm-ion 8
      :internal-ion 9
      :internal-17 20
      :internal-18 11
      }
     )
