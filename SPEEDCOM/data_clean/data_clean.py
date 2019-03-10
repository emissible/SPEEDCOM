"""
Functions for assisting in the cleaning of the obtained information.
"""

def remove_deliminators(my_strings):
  """
  Remove deliminators from numbers (comma) so as
  to be able to process numbers as int or float types in place of
  strings.  Takes in 'my_array', an array of strings, and will return
  an array of floats.
  """
  my_array = []
  for i in my_strings:
    if ',' in i:
      tmp = i.split(",")
      i = tmp[0] + tmp[1]
    try:
      my_array.append(float(i))  
    except:
      print("String" + i + " not able to be cast to float, characters \
      other than ',' or '.'?")
  return np.asarray(my_array)
