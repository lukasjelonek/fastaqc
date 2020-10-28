from enum import Enum

# source https://www.ddbj.nig.ac.jp/ddbj/code-e.html
_dna_unamb_chars = ['A','T','G','C']
_rna_unamb_chars = ['A','U','G','C']
_aa_unamb_chars = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
_aa_nc_unamb_chars = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','O','U']
_dna_amb_chars = ['N','R','Y','S','W','K','M','B','D','H','V']
_rna_amb_chars = _dna_unamb_chars
_aa_amb_chars = ['X','B','Z','J']
_dna_special_chars = ['-'] # gap
_rna_special_chars = _dna_special_chars
_aa_special_chars = ['*', '-'] # stop and gap

class Alphabet(Enum):
  DNA = (_dna_unamb_chars, _dna_amb_chars, _dna_special_chars)
  RNA = (_rna_unamb_chars, _rna_amb_chars, _rna_special_chars)
  AA = (_aa_unamb_chars, _aa_amb_chars, _aa_special_chars )
  AA_NC = (_aa_nc_unamb_chars, _aa_amb_chars, _aa_special_chars)
  # other contains all characters from above, but everything from the ambiguous section
  # is added to the non-ambiguous section
  OTHER = (
    list(set(_dna_unamb_chars + _rna_unamb_chars + _aa_unamb_chars + _dna_amb_chars + _rna_amb_chars )),
    _aa_amb_chars,
    list(set(_dna_special_chars + _rna_special_chars + _aa_special_chars))
    )

  def __init__(self, unambiguous_chars, ambiguous_chars, special_chars):
    self.unambiguous_chars = unambiguous_chars
    self.ambiguous_chars = ambiguous_chars
    self.special_chars = special_chars
    self.all_chars = unambiguous_chars + ambiguous_chars + special_chars
