
def mendels_probability(k: int, m: int, n: int) -> float:
    total = k + m + n    
    k_n = (k / total) * (n / (total - 1)) + (n / total) * (k / (total - 1))
    k_k = (k / total) * ((k - 1) / (total - 1))
    k_m = (k / total) * (m / (total - 1)) + (m / total) * (k / (total - 1))
    m_m = (m / total) * ((m - 1) / (total - 1))
    m_n = (m / total) * (n / (total - 1)) + (n / total) * (m / (total - 1))
    return k_k + k_n + k_m + (0.5 * m_n) + (0.75 * m_m)