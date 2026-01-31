# circulant_core_hadamard_search.py
#
# We search for a binary "core row" c of length m = 4n-1 with exactly 2n ones
# (equivalently: 2n entries equal to -1 under (-1)^bit, and 2n-1 entries +1).
#
# A normalized Hadamard with *circulant core* corresponds to:
# for every nonzero shift s (mod m), the dot product of c with shift(c,s) is -1.
# In 0/1 language this is equivalent to:
#   #agreements(c, shift(c,s)) == 2n-1
# for every s = 1..m-1.

def agreements_with_shift(bits, s):
    """Number of indices j with bits[j] == bits[(j+s) % m]."""
    m = len(bits)
    return sum(1 for j in range(m) if bits[j] == bits[(j + s) % m])


def is_circulant_core_candidate(bits, n):
    m = len(bits)                 # should be 4n-1
    target = 2 * n - 1
    # shift 0 is trivial; check all nonzero shifts modulo m
    for s in range(1, m):
        if agreements_with_shift(bits, s) != target:
            return False
    return True


def next_fixed_weight_lex(bits):
    """
    Advance to the next lexicographic (0<1) bitstring with the same Hamming weight.
    Returns True if advanced, False if already at the last one.
    """
    m = len(bits)

    # Find rightmost occurrence of '01'
    # (i.e., bits[i]==0 and bits[i+1]==1). Flip it to '10' and pack remaining 1s to the right.
    for i in range(m - 2, -1, -1):
        if bits[i] == 0 and bits[i + 1] == 1:
            ones_to_right = sum(bits[i + 1:])  # includes the 1 at i+1
            bits[i] = 1
            # clear suffix
            for k in range(i + 1, m):
                bits[k] = 0
            # place remaining ones (ones_to_right - 1) at the far right
            for k in range(m - (ones_to_right - 1), m):
                bits[k] = 1
            return True

    return False  # no '01' found => we're done


def search(n):
    m = 4 * n - 1
    weight = 2 * n

    # Lexicographically first bitstring of length m with exactly 'weight' ones:
    # (m-weight) zeros followed by 'weight' ones.
    bits = [0] * (m - weight) + [1] * weight

    while True:
        if is_circulant_core_candidate(bits, n):
            return bits  # found

        if not next_fixed_weight_lex(bits):
            return None  # exhausted


if __name__ == "__main__":
    n = int(input("n = ").strip())
    sol = search(n)

    if sol is None:
        print(f"No circulant core Hadamard found (by this brute-force search) for order {4*n}.")
    else:
        row = "".join(map(str, sol))
        print(f"Found circulant core row for order {4*n}: {row}")
        # If you want the ±1 version using (-1)^bit:
        # print("±1 core:", " ".join(str(1 if b == 0 else -1) for b in sol))