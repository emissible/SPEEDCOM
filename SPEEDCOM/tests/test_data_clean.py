"""
Unit tests for the data_clean.py file
"""

import unittest

from SPEEDCOM.data_clean import data_clean 


try:
  NotADirectoryError
except NameError:
  NotADirectoryError = OSError

# Unit tests                                                                    
class remove_deliminators(unittest.TestCase):
  # Assert that number without deliminators returns same number
  def test_nondeliminated_number(self):
    """
    Tests to ensure a number without a separator will return.
    """
    assert remove_deliminators(['42'])[0] == 42, \
      "remove_deliminators unable to handle non-deliminated numbers"
    #return

  def test_comma_delim_number(self):
    """
    Tests the removal of a comma deliminator.
    """
    assert remove_deliminators(['4,242'])[0] == 4242, \
      "remove_deliminators unable to remove comma deliminator"
    return

  def test_decimal_delim_number(self):
    """
    Tests a decimal is still able to return (python cast to float should
    be able to handle this by default.
    """
    assert remove_deliminators(['42.42'])[0] == 42.42, \
      'remove_deliminators unable to handle decimal point in number'
    return

  def test_decimal_and_comma(self):
    """
    Tests a number deliminated with a comma and containing a decimal will
    be able to be properly seperated and returned.
    """
    assert remove_deliminators(['4,242.42'])[0] == 4242.42, \
      'remove_deliminators unable to handle commas and decimals'
    return
  
  def test_nonnumber(self):
    """
    Tests that a non-numerical string does not return a number.  Should be
    handled by the type casting in python.
    """
    self.assertRaises(Exception, \
      lambda:remove_deliminators('NOT_A_NUMBER,But.has COMMA'))
    return
