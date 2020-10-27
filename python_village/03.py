"""
Problem
Given: A string s of length at most 200 letters and four integers a, b, c and d.

Return: The slice of this string from indices a through b and c through d (with space in between), inclusively. In other words, we should include elements s[b] and s[d] in our slice.
"""

s = "mxk2cEEcvLoY8b8TzMqW6iFQRUmxWk1BDS20xHAQizwCorytophanesxi1blYtycirciaRlqMrZSSlRuE7WTJkS9wqTc50uVsau0E1gnysqbKljPsSjk4VRw1fb1pteqvBKrU9jLMA0UMgXlnRjEhFcbnsjyvWXSz73D."

a = 43
b = 54
c = 63
d = 68


first = s[a:b + 1]

second = s[c:d + 1]

print(first + " " + second)

# Corytophanes circia