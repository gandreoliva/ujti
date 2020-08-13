import re
import sys
from textwrap import wrap

"""
Auxiliary script: changes single precision to double precision

usage: python f90-single2double.py file.f90

Intended use: The output from Maxima sometimes (Lisp dependent) contains
single precision expressions instead of double precision. For example,
    1.0/2.0
but for its use in Fortran, the form
    1d0/2d0
is necessary. This script does the changes.
"""


file = open(sys.argv[1],'r')
text = file.read()
file.close()

file = open(sys.argv[1],'w')

# Gets rid of the continuation lines
text = re.sub(r"&\s*\n\s*&","",text)

# Changes decimal numbers to real8 notation
# case 1.03 -> 1.03d0
text = re.sub(r"\.(?=(\d+))\1(?!e|d)",r".\1e0",text)
r"""
Explanation:
\. : match the dot
(?=...) : folllowed by
(?=(...))\1 : workaround to atomic groups (not supported by re)
    looks for a whole match ("give me all or nothing")
(?!e|d) : not followed by e (or a "d", if already in double form)
(?=(\d+))\1(?!e) : workaround, since python doesn't support
    possesive quantifiers. Matches a whole seq. of digits \d+ that
    are not followe by e. Example:
    (\d+)(?!e) on 1.832e92+9.03 matches 1.83 and 9.03 but
    (?=(\d+))\1(?!e) does not match the first number and only
    matches 9.03
"""

# Changes single precision exponential notation to double precision exp. notation
# case 2.6e9 -> 2.6d9
text = re.sub(r"\.(\d+)e",r".\1d", text)
r"""
Explanation:
\. : matches the dot
(\d+) : matches a sequence of digits
e : followed by e
"""



# Wraps the lines with Fortran continuation lines
# text = "&\n&".join(wrap(text,72))
text = re.sub(r"(.{72})",r"\1&\n&",text)


file.write(text)
file.close()
