import argparse
import sys

ALPHABETS = {
    "ru": "абвгдежзийклмнопрстуфхцчшщъыьэюя",
    "en": "abcdefghijklmnopqrstuvwxyz", 
}


def load_text(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def save_text(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def load_key_letter(path: str, alphabet: str) -> str:
    raw = load_text(path).strip()
    if not raw:
        raise ValueError("Файл ключа пуст: нужна одна начальная буква K0.")
    k0 = raw[0].lower()
    if k0 not in alphabet:
        raise ValueError(f"Начальная буква ключа '{k0}' не входит в выбранный алфавит.")
    return k0


def vigenere_autokey_encrypt(text: str, k0: str, alphabet: str) -> str:
    m = len(alphabet)
    idx = {ch: i for i, ch in enumerate(alphabet)}
    prev_key_idx = idx[k0]
    out = []
    for ch in text:
        low = ch.lower()
        if low in idx:
            p = idx[low]
            c = (p + prev_key_idx) % m
            out_ch = alphabet[c]
            if ch.isupper():
                out_ch = out_ch.upper()
            out.append(out_ch)
            prev_key_idx = c
        else:
            out.append(ch)
    return "".join(out)


def vigenere_autokey_decrypt(text: str, k0: str, alphabet: str) -> str:
    m = len(alphabet)
    idx = {ch: i for i, ch in enumerate(alphabet)}
    prev_key_idx = idx[k0]
    out = []
    for ch in text:
        low = ch.lower()
        if low in idx:
            c = idx[low]
            p = (c - prev_key_idx) % m
            out_ch = alphabet[p]
            if ch.isupper():
                out_ch = out_ch.upper()
            out.append(out_ch)
            prev_key_idx = c
        else:
            out.append(ch)
    return "".join(out)


def main():
    parser = argparse.ArgumentParser(description="Шифр Вижинера с самоключом по шифротексту")
    parser.add_argument("mode", choices=["encrypt", "decrypt"])
    parser.add_argument("--in", dest="infile", required=True, help="входной текстовый файл")
    parser.add_argument("--key", dest="keyfile", required=True, help="файл с начальной буквой ключа K0")
    parser.add_argument("--out", dest="outfile", required=True, help="выходной текстовый файл")
    parser.add_argument("--alphabet", choices=["ru", "en"], default="ru", help="алфавит сообщения")
    args = parser.parse_args()

    alphabet = ALPHABETS[args.alphabet]
    k0 = load_key_letter(args.keyfile, alphabet)
    text = load_text(args.infile)

    if args.mode == "encrypt":
        result = vigenere_autokey_encrypt(text, k0, alphabet)
    else:
        result = vigenere_autokey_decrypt(text, k0, alphabet)

    save_text(args.outfile, result)
    print(f"Готово: {args.outfile}")


if __name__ == "__main__":
    main()
