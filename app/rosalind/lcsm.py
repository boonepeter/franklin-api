

from typing import List


def substr_in_all(arr: List[str], part: str) -> bool:
  for dna in arr:
    if dna.find(part) == -1:
      return False
  return True

def common_substr(arr: List[str], l: int) -> str:
  first = arr[0]
  for i in range(len(first) - l + 1):
    part = first[i: i + l]
    if substr_in_all(arr, part):
      return part
  return ""

def longest_common_substr(arr: List[str]) -> str:
  l = 0
  r = len(arr[0])
  while l + 1 < r:
    mid = (l + r) // 2
    if common_substr(arr, mid) != "":
      l = mid
    else:
      r = mid

  return common_substr(arr, l)