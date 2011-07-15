(ns glutton.mass
  "Important physical constants")

(def amino-acid-dictionary
     {:. {:average-mass 0, :monoisotopic-mass 0}
      :A {:average-mass 71.0788, :monoisotopic-mass 71.03711}
      :C {:average-mass 103.1448, :monoisotopic-mass 103.00919}
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
      :Y {:average-mass 163.176, :monoisotopic-mass 163.06333}})

(def CHON
  "Carbon, Hydrogen, Oxygen, and Nitrogen ([CHON][1]) mass constants.  These are important
  for computing masses of organic compounds.

  [1]: http://en.wikipedia.org/wiki/CHON "
     {:carbon {:monoisotopic-mass 12.00000000
               :average-mass 12.01073590}
      :hydrogen {:monoisotopic-mass 1.00782504
                 :average-mass 1.00794076}
      :oxygen {:monoisotopic-mass 15.99491463
               :average-mass 15.99940494}
      :nitrogen {:monoisotopic-mass 14.00307400
                 :average-mass 14.00674309}})

(def water-mass
  "Mass calculations for water (H<sub>2</sub>0)."
  {:monoisotopic-mass (+ (:monoisotopic-mass (:oxygen CHON))
                         (* 2 (:monoisotopic-mass (:hydrogen CHON))))
   :average-mass (+ (:average-mass (:oxygen CHON))
                    (* 2 (:average-mass (:hydrogen CHON))))})

(defn aa-mass
  "Retrieve the mass of a given `amino-acid`"
  [amino-acid mass-type]
  (get-in amino-acid-dictionary [amino-acid mass-type]))
