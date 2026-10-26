import math
from scipy.special import gammaincc

class LFSR:
    def __init__(self, n, taps, seed):
        self.n = n
        self.taps = taps  # множество индексов j с a_j = 1
        self.state = list(seed)
        assert len(self.state) == n
        assert any(self.state), "начальное состояние не должно быть нулевым"

    def step(self):
        out = self.state[0]
        fb = 0
        for j in self.taps:
            fb ^= self.state[j]
        self.state = self.state[1:] + [fb]
        return out

    def generate(self, length):
        return ''.join(str(self.step()) for _ in range(length))

    def period(self):
        start = tuple(self.state)
        seen = {start: 0}
        t = 0
        while True:
            self.step()
            t += 1
            cur = tuple(self.state)
            if cur in seen:
                return t - seen[cur]
            seen[cur] = t


def taps_from_poly(exponents, n):
    return sorted(e for e in exponents if e < n)



def check_polynomials():
    n = 6
    polys = {
        "f1 = x^6+x^3+1":   {6, 3, 0},
        "f2 = x^6+x^2+x+1": {6, 2, 1, 0},
        "f3 = x^6+x+1":     {6, 1, 0},
    }
    print("=== Проверка периода LFSR для трёх многочленов степени n=6 ===")
    print(f"Максимально возможный период T_max = 2^{n}-1 = {2**n - 1}\n")
    results = {}
    for name, exps in polys.items():
        taps = taps_from_poly(exps, n)
        seed = [0] * (n - 1) + [1]
        lf = LFSR(n, taps, seed)
        T = lf.period()
        is_primitive = (T == 2**n - 1)
        results[name] = (taps, T, is_primitive)
        print(f"{name:22s} taps(j<{n})={taps}  период T={T:3d}  "
              f"{'ПРИМИТИВЕН' if is_primitive else 'НЕ примитивен'}")
    print()
    return results


def longest_run_of_ones_test(bits):
    """Тест 'Самая длинная серия единиц в блоке', параметры для n>=6272."""
    n = len(bits)
    assert n >= 128
    if n < 6272:
        M, K, N = 8, 3, n // 8
        pi = [0.2148, 0.3672, 0.2305, 0.1875]
        vmap = lambda r: 0 if r <= 1 else (1 if r == 2 else (2 if r == 3 else 3))
    elif n < 750000:
        M, K, N = 128, 5, n // 128
        pi = [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124]
        vmap = lambda r: 0 if r <= 4 else (1 if r == 5 else (2 if r == 6 else (3 if r == 7 else (4 if r == 8 else 5))))
    else:
        M, K, N = 10000, 6, n // 10000
        pi = [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]
        vmap = lambda r: (0 if r <= 10 else (1 if r == 11 else (2 if r == 12 else
                          (3 if r == 13 else (4 if r == 14 else (5 if r == 15 else 6))))))

    v = [0] * (K + 1)
    for i in range(N):
        block = bits[i * M:(i + 1) * M]
        maxrun = run = 0
        for ch in block:
            if ch == '1':
                run += 1
                maxrun = max(maxrun, run)
            else:
                run = 0
        v[vmap(maxrun)] += 1

    chi2 = sum((v[i] - N * pi[i]) ** 2 / (N * pi[i]) for i in range(K + 1))
    p = gammaincc(K / 2, chi2 / 2)
    return p, chi2, v


def non_overlapping_template_matching_test(bits, template="000000001", M=125000):
    """Тест на совпадение непересекающихся шаблонов."""
    n = len(bits)
    m = len(template)
    N = n // M
    mu = (M - m + 1) / 2 ** m
    var = M * (1.0 / 2 ** m - (2 * m - 1) / 2 ** (2 * m))

    W = []
    for i in range(N):
        block = bits[i * M:(i + 1) * M]
        count = 0
        pos = 0
        while pos <= M - m:
            if block[pos:pos + m] == template:
                count += 1
                pos += m
            else:
                pos += 1
        W.append(count)

    chi2 = sum((w - mu) ** 2 / var for w in W)
    p = gammaincc(N / 2, chi2 / 2)
    return p, chi2, W


def berlekamp_massey(bits):
    n = len(bits)
    s = [int(b) for b in bits]
    c = [0] * n
    b = [0] * n
    c[0] = 1
    b[0] = 1
    L = 0
    m = -1
    for N in range(n):
        d = s[N]
        for i in range(1, L + 1):
            d ^= c[i] & s[N - i]
        if d == 1:
            t = c.copy()
            shift = N - m
            for i in range(n - shift):
                c[i + shift] ^= b[i]
            if L <= N / 2:
                L = N + 1 - L
                m = N
                b = t
    return L


def linear_complexity_test(bits, M=500):
    n = len(bits)
    N = n // M
    pi = [0.010417, 0.031250, 0.125000, 0.500000, 0.250000, 0.062500, 0.020833]

    mu = M / 2 + (9 + (-1) ** (M + 1)) / 36 - (M / 3 + 2 / 9) / 2 ** M

    v = [0] * 7
    Ls = []
    for i in range(N):
        block = bits[i * M:(i + 1) * M]
        L = berlekamp_massey(block)
        Ls.append(L)
        T = ((-1) ** M) * (L - mu) + 2 / 9
        if T <= -2.5:
            v[0] += 1
        elif T <= -1.5:
            v[1] += 1
        elif T <= -0.5:
            v[2] += 1
        elif T <= 0.5:
            v[3] += 1
        elif T <= 1.5:
            v[4] += 1
        elif T <= 2.5:
            v[5] += 1
        else:
            v[6] += 1

    chi2 = sum((v[i] - N * pi[i]) ** 2 / (N * pi[i]) for i in range(7))
    p = gammaincc(3, chi2 / 2)
    return p, chi2, v


def load_bits_from_file(path, length=1_000_000):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        raw = f.read()
    digits = ''.join(ch for ch in raw if ch in '01')
    return digits[:length]


if __name__ == "__main__":
    results = check_polynomials()
    primitive = [name for name, (taps, T, ok) in results.items() if ok]
    print("Примитивные многочлены по модулю 2:", primitive, "\n")

    # --- Шаг 3a: валидация тестов на бинарных разложениях числа pi и e ---
    pi_bits = load_bits_from_file("/mnt/user-data/uploads/data.pi")
    e_bits = load_bits_from_file("/mnt/user-data/uploads/data.e")
    print(f"Загружено бит: pi={len(pi_bits)}, e={len(e_bits)}\n")

    print("=== Валидация тестов (ожидаемые p-value из таблицы методички) ===")
    for name, bits in [("pi", pi_bits), ("e", e_bits)]:
        p1, chi1, _ = non_overlapping_template_matching_test(bits)
        p2, chi2, _ = linear_complexity_test(bits)
        p3, chi3, _ = longest_run_of_ones_test(bits)
        print(f"[{name}] Непересек.шаблонов p={p1:.6f} | "
              f"Лин.сложность p={p2:.6f} | Самая длинная серия единиц p={p3:.6f}")
    print("Ожидалось (pi): 0.078790 | 0.826335 | 0.718945")
    print("Ожидалось (e):  0.165757 | 0.255475 | 0.024390\n")

    n = 6
    taps = taps_from_poly({6, 3, 0}, n)
    seed = [0, 0, 0, 0, 0, 1]
    lf = LFSR(n, taps, seed)
    seq = lf.generate(1_000_000)
    print(f"Сгенерирована LFSR-последовательность, длина={len(seq)}, период={2**n-1} "
          f"(повторяется {1_000_000 // (2**n - 1)} раз, т.к. n мало)\n")

    alpha = 0.01
    print("=== Результаты тестирования сгенерированной LFSR-последовательности (alpha=0.01) ===")
    p1, chi1, _ = non_overlapping_template_matching_test(seq)
    p2, chi2, _ = linear_complexity_test(seq)
    p3, chi3, _ = longest_run_of_ones_test(seq)
    for label, p in [("Непересекающихся шаблонов", p1),
                      ("На линейную сложность", p2),
                      ("Самая длинная серия единиц в блоке", p3)]:
        verdict = "СЛУЧАЙНА (тест пройден)" if p >= alpha else "НЕ случайна (тест НЕ пройден)"
        print(f"  {label:38s} p-value={p:.6f}  -> {verdict}")
