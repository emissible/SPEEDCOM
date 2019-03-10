"""
Unit tests for the data_clean.py file
"""

import unittest

import ../data_clean/data_clean


# Unit tests                                                                    
class remove_deliminators(unittest.TestCase):
  # Assert that number without deliminators returns same number
  def test_nondeliminated_number(self):
    """
    Tests to ensure a number without will return.
    """
    assert int(data_clean.remove_deliminators(['42']) == 42, \
      'remove_deliminators unable to handle non-deliminated numbers'
    return

  def test_comma_delim_number(unittest.TestCase):
    """
    Tests the removal of a comma deliminator.
    """
    assert int(data_clean.remove_deliminators(['4,242'])) == 4,242, \
      'remove_deliminators unable to remove comma deliminator'
    return

  def test_decimal_delim_number(unittest.TestCase):
    """
    Tests a decimal is still able to return (python cast to float should
    be able to handle this by default.
    """
    assert (data_clean.remove_deliminators(['42.42']) == 42.42, \
      'remove_deliminators unable to handle decimal point in number'
    return

  def test_decimal_and_comma(unittest.TestCase):
    """
    Tests a number deliminated with a comma and containing a decimal will
    be able to be properly seperated and returned.
    """
    assert data_clean.remove_deliminators(['4,242.42']) == 4242.42, \
      'remove_deliminators unable to handle commas and decimals'
    return
  
  def test_nonnumber(self):
    """
    Tests that a non-numerical string does not return a number.  Should be
    handled by the type casting in python.
    """
    self.assertRaises(Exception, \
      lambda:data_clean.remove_deliminators('NOT_A_NUMBER,But.has COMMA'))
    return
