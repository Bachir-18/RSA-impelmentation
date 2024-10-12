#include "Bignum.hpp"
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

Bignum::Bignum(int x) {
  if (x >= 0) {
    tab = vector<uint32_t>(1, x);
    isPositive = true;
  } else {
    tab = vector<uint32_t>(1, -x);
    isPositive = false;
  }
}

Bignum::Bignum(unsigned x) {
  tab = vector<uint32_t>(1, x);
  isPositive = true;
}

Bignum::Bignum(string s) {
  // If s = '' return Bignum(0)
  if (s.size() == 0) {
    tab = vector<uint32_t>(1, 0);
    isPositive = true;
  } else {
    // Handle negative numbers and delete the '-' sign
    if (s.at(0) == '-') {
      s.erase(0, 1);
      isPositive = false;
    } else
      isPositive = true;
    tab = vector<uint32_t>(1, 0);
    unsigned pos = 0;   // Store the leftmost bit position
    unsigned block = 0; // Store the leftmost block of the tab vector
    // delete unsignificant zeros if there are any
    while (s.size() != 0 && s.at(0) == '0')
      s.erase(0, 1);
    // Compute the number in base pow(2,32)
    while (s.size() != 0) {
      string tmp = "";
      unsigned n = 0;
      for (auto c : s) {
        n *= 10;
        n += c - '0';
        tmp += to_string(n / 2);
        n = n % 2;
      }
      s = move(tmp);
      while (s.size() != 0 && s.at(0) == '0')
        s.erase(0, 1);
      tab[block] |= n << pos;
      // If the block contain less than 32 bits increment pos
      if (pos < 31)
        pos += 1;
      // Otherwise set pos to zero and go to another block
      else {
        pos = 0;
        block += 1;
        tab.emplace_back(0);
      }
    }
  }
}

// Print Bignum in hexadecimal base
void Bignum::printHex() const {
  if (!isPositive)
    cout << "-";
  for (auto rit = tab.rbegin(); rit != tab.rend(); ++rit) {
    cout << hex;
    if (rit != tab.rbegin())
      cout << setfill('0') << setw(8);
    cout << *rit << dec;
  }
}

// Function to add two big integers with the same sign
Bignum addSameSign(Bignum const &x, Bignum const &y) {
  if (y.tab.size() > x.tab.size())
    return addSameSign(y, x); // Swap x and y if y has more digits

  Bignum z(0u);
  z.tab.resize(x.tab.size());
  uint64_t c = 0; // Carry

  // Add corresponding digits of x and y
  for (unsigned i = 0; i < y.tab.size(); ++i) {
    uint64_t tmp = x.tab[i] + (y.tab[i] + c);
    z.tab[i] = tmp;
    c = tmp >> 32; // Store the carry for the next digit addition
  }

  // Add the remaining digits of x
  for (unsigned i = y.tab.size(); i < x.tab.size(); ++i) {
    uint64_t tmp = x.tab[i] + c;
    z.tab[i] = tmp;
    c = tmp >> 32; // Store the carry for the next digit addition
  }

  // If there is a carry after all additions, add it as a new digit
  if (c != 0)
    z.tab.emplace_back(1);

  return z;
}

// Function to calculate the greatest common divisor (gcd) of two big integers
Bignum pgcd(const Bignum &x, const Bignum &y) {
  Bignum zero = Bignum(0);
  Bignum u = x;
  Bignum v = y;
  Bignum tmp;

  // Euclidean algorithm to find the gcd
  while (v != zero) {
    tmp = v;
    v = u % v;
    u = tmp;
  }

  return u;
}

// Function to calculate the modular inverse of a big integer
Bignum inverseMod(const Bignum &x, const Bignum &y) {
  // We compute the inverse modulo by using the Extended Euclidean Algorithm
  // We assume that y is a prime number and x is greater than zero
  Bignum x0 = Bignum(1);
  Bignum y0 = Bignum(0);
  Bignum x1 = Bignum(0);
  Bignum y1 = Bignum(1);
  Bignum tmp_x;
  Bignum tmp_y;
  Bignum q0;
  Bignum r0;
  Bignum gcd = x, b0 = y;

  while (b0 != 0) {
    q0 = gcd / b0;
    r0 = gcd % b0;
    gcd = b0;
    b0 = r0;

    tmp_x = x0;
    x0 = x1;
    x1 = tmp_x - (q0 * x1);

    tmp_y = y0;
    y0 = y1;
    y1 = tmp_y - (q0 * y1);
  }

  // If the gcd is not equal to 1, the modulus value is not a prime number or x
  // is equal to zero
  if (gcd != 1) {
    cerr << "No inverse found for this chosen number" << endl;
    exit(-1);
  }
  return x0 % y;
}

// Function to subtract two big integers (assume x >= y >= 0)
Bignum SubtractX_Y(Bignum const &x, Bignum const &y) {
  Bignum z(0u);
  z.tab.resize(x.tab.size());
  uint64_t c = 0; // Borrow

  // Subtract corresponding digits of x and y
  for (unsigned i = 0; i < y.tab.size(); ++i) {
    uint64_t tmp = x.tab[i] - (y.tab[i] + c);
    z.tab[i] = tmp;
    if ((tmp >> 32) != 0)
      c = 1; // Set the borrow for the next digit subtraction
    else
      c = 0;
  }

  // Subtract the remaining digits of x
  for (unsigned i = y.tab.size(); i < x.tab.size(); ++i) {
    uint64_t tmp = x.tab[i] - c;
    z.tab[i] = tmp;
    if ((tmp >> 32) != 0)
      c = 1; // Set the borrow for the next digit subtraction
    else
      c = 0;
  }

  return z;
}

// Function to compare the absolute values of two big integers
bool compareAbs(Bignum const &x, Bignum const &y) {
  unsigned nx = x.tab.size();
  while (nx > 0 && x.tab[nx - 1] == 0)
    --nx;
  unsigned ny = y.tab.size();
  while (ny > 0 && y.tab[ny - 1] == 0)
    --ny;

  if (nx != ny)
    return nx > ny; // Compare the number of digits

  // Compare each digit from most significant to least significant
  while (nx > 0) {
    --nx;
    if (x.tab[nx] != y.tab[nx])
      return x.tab[nx] > y.tab[nx]; // Compare the digits
  }

  return true; // The absolute values are equal
}
// Addition operator overload
Bignum operator+(Bignum const &x, Bignum const &y) {
  if (x.isPositive == y.isPositive) {
    auto z = addSameSign(x, y);
    z.isPositive = x.isPositive;
    return z;
  } else {
    if (x.isPositive) {
      if (compareAbs(x, y))
        return SubtractX_Y(x, y);
      auto z = SubtractX_Y(y, x);
      z.isPositive = false;
      return z;
    } else {
      if (compareAbs(y, x))
        return SubtractX_Y(y, x);
      auto z = SubtractX_Y(x, y);
      z.isPositive = false;
      return z;
    }
  }
}

// Subtract operator overload
Bignum operator-(Bignum const &x, Bignum const &y) {
  if (x.isPositive == y.isPositive) {
    if (compareAbs(x, y)) {
      auto z = SubtractX_Y(x, y);
      z.isPositive = x.isPositive;
      return z;
    } else {
      auto z = SubtractX_Y(y, x);
      z.isPositive = !x.isPositive;
      return z;
    }
  } else {
    auto z = addSameSign(x, y);
    z.isPositive = x.isPositive;
    return z;
  }
}

// Multiply operator overload
Bignum operator*(Bignum const &x, Bignum const &y) {
  Bignum z;
  int t = x.tab.size() + y.tab.size();
  int n = x.tab.size();
  int m = y.tab.size();
  uint64_t b = static_cast<uint64_t>(1) << 32;
  uint64_t c;
  uint64_t temp;

  // Initialize z with zeros
  for (int i = 0; i < t; i++) {
    z.tab.push_back(0);
  }

  // Perform the multiplication
  for (int i = 0; i < n; i++) {
    c = 0;
    for (int j = 0; j < m; j++) {
      temp = (uint64_t)z.tab[i + j] +
             ((uint64_t)x.tab[i] * (uint64_t)y.tab[j]) + c;
      z.tab[i + j] = temp % b;
      c = temp / b;
    }
    z.tab[i + m] += c;
  }

  // Set the sign of the result
  z.isPositive = !(x.isPositive ^ y.isPositive);

  // Remove leading zeros from the result
  int j = z.tab.size() - 1;
  while (j > 0 && z.tab[j] == 0) {
    z.tab.pop_back();
    j--;
  }

  return z;
}

// Division function
pair<Bignum, Bignum> division(Bignum const &x, Bignum const &y) {
  if (y == Bignum(0)) {
    cerr << "Error: Division by 0!" << endl;
    exit(1);
  }

  Bignum a = x, b = y;
  Bignum r = x;
  Bignum q = Bignum(0);
  Bignum tmp;

  // Handle special cases
  if (a < b) {
    return make_pair(Bignum(0), Bignum(x));
  }
  if (a == b) {
    return make_pair(Bignum(1), Bignum(0));
  }

  // Normalize the divisor
  tmp = b;
  while (tmp < a) {
    tmp = tmp << 1;
  }
  if (tmp > a) {
    tmp = tmp >> 1;
  }

  // Perform the division
  while (tmp >= b) {
    q = q << 1;

    if (r >= tmp) {
      q = q + 1;
      r = r - tmp;
    }

    tmp = tmp >> 1;
  }

  return make_pair(q, r);
}

// Left bit shift operator
Bignum &Bignum::operator<<=(unsigned n) {
  if (n == 0)
    return *this;
  unsigned n_bits = n % 32;
  unsigned n_blocks = n / 32;

  // Shift the digits and adjust the size
  if (n_bits == 0) {
    tab.resize(tab.size() + n_blocks);
    for (unsigned i = tab.size(); i-- != n_blocks;)
      tab[i] = tab[i - n_blocks];
  } else {
    tab.resize(tab.size() + n_blocks + 1);
    for (unsigned i = tab.size(); i-- != n_blocks + 1;)
      tab[i] = (tab[i - n_blocks] << n_bits) |
               (tab[i - (n_blocks + 1)] >> (32 - n_bits));
    tab[n_blocks] = tab[0] << n_bits;
  }

  // Remove leading zeros
  unsigned ntab = tab.size();
  while (ntab > 1 && tab[ntab - 1] == 0)
    --ntab;
  tab.resize(ntab);

  return *this;
}

// Right bit shift operator
Bignum &Bignum::operator>>=(unsigned n) {
  if (n == 0)
    return *this;
  unsigned n_bits = n % 32;
  unsigned n_blocks = n / 32;

  // Shift the digits and adjust the size
  if (n_bits == 0) {
    for (unsigned i = 0; i + n_blocks < tab.size(); ++i)
      tab[i] = tab[i + n_blocks];
    if (tab.size() > n_blocks)
      tab.resize(tab.size() - n_blocks);
    else {
      tab.resize(1);
      tab[0] = 0;
    }
  } else {

    for (unsigned i = 0; i + 1 + n_blocks < tab.size(); ++i)
      tab[i] = (tab[i + n_blocks] >> n_bits) |
               (tab[i + (n_blocks + 1)] << (32 - n_bits));
    tab[tab.size() - n_blocks - 1] = tab[tab.size() - 1] >> n_bits;
  }

  // Remove leading zeros
  unsigned ntab = tab.size();
  while (ntab > 1 && tab[ntab - 1] == 0)
    --ntab;
  tab.resize(ntab);

  return *this;
}
// Comparison: Less than operator
bool operator<(Bignum const &x, Bignum const &y) {
  if (x.isPositive != y.isPositive)
    return y.isPositive; // If signs are different, return the comparison based
                         // on y's sign
  unsigned n = x.size() - 1, t = y.size() - 1;
  while (n > 0 && x[n] == 0)
    --n; // Find the highest non-zero element index in x
  while (t > 0 && y[t] == 0)
    --t; // Find the highest non-zero element index in y
  if (n != t)
    return x.isPositive ^ (n > t); // If sizes are different, return the
                                   // comparison based on sizes and signs
  while (n != 0) {
    if (x[n] != y[n])
      return x.isPositive ^
             (x[n] > y[n]); // Compare non-zero elements from highest to lowest
    --n;
  }
  if (x[0] == y[0])
    return false; // If all elements are equal, return false
  return x.isPositive ^
         (x[0] > y[0]); // Compare the lowest element and return based on sign
}

// Comparison: Less than or equal to operator
bool operator<=(Bignum const &x, Bignum const &y) {
  if (x.isPositive != y.isPositive)
    return y.isPositive; // If signs are different, return the comparison based
                         // on y's sign
  unsigned n = x.size() - 1, t = y.size() - 1;
  while (n > 0 && x[n] == 0)
    --n; // Find the highest non-zero element index in x
  while (t > 0 && y[t] == 0)
    --t; // Find the highest non-zero element index in y
  if (n != t)
    return x.isPositive ^ (n > t); // If sizes are different, return the
                                   // comparison based on sizes and signs
  while (n != 0) {
    if (x[n] != y[n])
      return x.isPositive ^
             (x[n] > y[n]); // Compare non-zero elements from highest to lowest
    --n;
  }
  if (x[0] == y[0])
    return true; // If all elements are equal, return true
  return x.isPositive ^
         (x[0] > y[0]); // Compare the lowest element and return based on sign
}

// Comparison: Equality operator
bool operator==(Bignum const &x, Bignum const &y) {
  unsigned n = x.size() - 1, t = y.size() - 1;
  while (n > 0 && x[n] == 0)
    --n; // Find the highest non-zero element index in x
  while (t > 0 && y[t] == 0)
    --t; // Find the highest non-zero element index in y
  if (n != t)
    return false; // If sizes are different, return false
  if (n != 0 && x.isPositive != y.isPositive)
    return false; // If signs are different and not zero, return false
  while (n != 0) {
    if (x[n] != y[n])
      return false; // Compare non-zero elements from highest to lowest
    --n;
  }
  return (x[0] == y[0]) && (x.isPositive == y.isPositive ||
                            x[0] == 0); // Compare the lowest element and signs
}

// Division operator
Bignum operator/(Bignum const &x, Bignum const &y) {
  auto xx = x;
  xx.isPositive = true;
  auto yy = y;
  yy.isPositive = true;
  auto p = division(move(xx), move(yy)); // Perform division
  if (x.isPositive != y.isPositive)
    p.first.isPositive = false; // Adjust sign based on x and y signs
  return p.first;
}

// Modulus operator
Bignum operator%(Bignum const &x, Bignum const &y) {
  auto xx = x;
  xx.isPositive = true;
  auto yy = y;
  yy.isPositive = true;
  auto p = division(move(xx), move(yy)); // Perform division
  if (!x.isPositive) {
    p.second.isPositive = false; // Adjust sign based on x sign
    if (y.isPositive)
      p.second += y; // Add y if y is positive
    else
      p.second -= y; // Subtract y if y is negative
  }
  return p.second;
}

// Resize the Bignum by removing leading zeros
void Bignum::resize() {
  // Remove leading zeros
  while (this->size() > 1 && this->tab.back() == 0) {
    this->tab.pop_back();
  }
  // If the Bignum is zero, set its sign to positive
  if (this->size() == 1 && this->tab[0] == 0) {
    this->isPositive = 1;
  }
}

// Output operator for Bignum
std::ostream &operator<<(std::ostream &flux, Bignum const &x) {
  if (x == Bignum(0))
    flux << 0; // If Bignum is zero, output 0
  else {
    if (!x.isPositive)
      flux << "-"; // Output '-' if Bignum is negative
    auto p = division(x, 10);
    p.first.isPositive = true;
    if (p.first != Bignum(0))
      flux << p.first;   // Output non-zero part
    flux << p.second[0]; // Output the lowest element
  }
  return flux;
}

// Modular exponentiation
Bignum modular_exponentiation(const Bignum &x, const Bignum &e,
                              const Bignum &m) {
  if (e == Bignum(0))
    return Bignum(1); // Base case: if exponent is 0, return 1

  Bignum res(1);
  Bignum a(x % m);
  Bignum exp(e);
  while (exp > 0) {
    if ((exp[0] & 1) == 1) {
      res = (res * a) % m; // Multiply and take modulus if exponent bit is 1
    }
    a = (a * a) % m; // Square and take modulus for the next iteration
    exp >>= 1;       // Right shift the exponent bits
  }

  return res;
}
// Generates a random odd Bignum of size n
Bignum odd_bignum_generator(unsigned int n) {
  unsigned int size = n / 32; // Number of 32-bit elements in the Bignum
  unsigned int bit = n % 32;  // Bit position within the last element
  Bignum x;
  uint64_t b = uint64_t(1) << 32;
  std::random_device rd; // Random device for seed generation
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint32_t> distribution(
      0, b - 1); // Uniform distribution

  if (bit == 0) {
    // If the bit position is 0, the last element is zero and the Bignum size is
    // 'size'
    x = Bignum(0);
    x.tab.resize(size);              // Resize the Bignum vector to 'size'
    unsigned nb = distribution(gen); // Generate a random number
    x[size - 1] = nb; // Assign the random number to the last element
    for (unsigned int i = 0; i < size - 1; i++) {
      x[i] = distribution(
          gen); // Generate random numbers for the remaining elements
    }
  } else {
    // If the bit position is non-zero, the Bignum size is 'size + 1'
    x = Bignum(0);
    x.tab.resize(size + 1); // Resize the Bignum vector to 'size + 1'
    x[size] =
        distribution(gen) &
        ((1 << bit) -
         1); // Generate a random number and mask it with the appropriate bits
    for (unsigned int i = 0; i < size; i++) {
      x[i] = distribution(
          gen); // Generate random numbers for the remaining elements
    }
  }

  if ((x % Bignum(2) == 0) && x != 2) {
    x = x + Bignum(1); // If x is even and not equal to 2, increment it by 1 to
                       // make it odd
  }

  return x; // Return the generated Bignum
}

// Performs the Fermat primality test on a given Bignum
bool FermatTest(Bignum n) {
  if (n == 2 || n == 3 || n == 5 || n == 7)
    return true; // n is a known prime number
  if (n == 1)
    return false;            // 1 is not a prime number
  int list[] = {2, 3, 5, 7}; // List of bases for the Fermat test
  for (int i = 0; i < 4; i++) {
    Bignum res = modular_exponentiation(Bignum(list[i]), n - 1,
                                        n); // Compute (list[i]^(n-1)) mod n
    if (res != Bignum(1)) {
      return false; // n is composite if the result is not 1
    }
  }
  return true; // n is likely a prime number
}

Bignum prime_generator(int n) {
  if (n <= Bignum(1)) {
    cerr << "The integer given must be bigger than 1" << endl;
    exit(1);
  }
  Bignum x = odd_bignum_generator(n);
  while (!FermatTest(x)) {
    x = x + Bignum(2);
  }
  return x;
}

vector<Bignum> key_generator(int x) {
  if (x < 10) {
    cerr << "Please choose big numbers in order to secure the encryption !"
         << endl;
    exit(1);
  }
  cout << "Generating keys of " << x << " bits" << endl;

  // Generate two prime numbers, p and q
  Bignum p = prime_generator(x);
  while (p == 1) {
    p = prime_generator(x);
  }
  Bignum q = prime_generator(x);
  while (q == p || q == 1) {
    q = prime_generator(x);
  }

  // Calculate n = p * q and phi = (p - 1) * (q - 1)
  Bignum n = p * q;
  Bignum phi = (p - Bignum(1)) * (q - Bignum(1));

  // Select an exponent e and calculate its modular inverse d
  Bignum e = Bignum(65537);
  Bignum d = inverseMod(e, phi);

  // Store the keys in a vector and return the result
  vector<Bignum> result = {e, d, n};
  return result;
}

vector<Bignum> encode(string text, int n) {
  // Encode the text using n-bit integers (Bignum)
  Bignum m(0);
  std::vector<Bignum> Vector;
  int j = 0;
  for (size_t i = 0; i < text.length(); i++) {
    m = m << 8;
    m = m + text[i];
    j = j + 1;
    if (j == n / 8) {
      Vector.push_back(m);
      m = 0;
      j = 0;
    }
  }
  if (m != 0) {
    Vector.push_back(m);
  }

  return Vector;
}

string decode(vector<Bignum> encodedText, int n) {
  // Decode the integers (Bignum) back to text
  string result;

  while (encodedText[0] > 0) {
    if (encodedText[encodedText.size() - 1] == 0) {
      encodedText.pop_back();
    }
    Bignum byte = encodedText[encodedText.size() - 1] % Bignum(256);
    result = static_cast<char>(byte.tab[0]) + result;
    encodedText[encodedText.size() - 1] =
        encodedText[encodedText.size() - 1] >> 8;
  }

  return result;
}

vector<Bignum> encrypt(vector<Bignum> &message, Bignum e, Bignum n) {
  // Encrypt the message using the exponent e and modulo n
  vector<Bignum> result;
  Bignum tmp(0);
  for (int i = 0; i < static_cast<int>(message.size()); i++) {
    tmp = modular_exponentiation(message[i], e, n);
    result.push_back(tmp);
  }
  return result;
}

vector<Bignum> decrypt(vector<Bignum> &ciphertext, Bignum d, Bignum n) {
  // Decrypt the ciphertext using the exponent d and modulo n
  vector<Bignum> result;
  Bignum tmp(0);
  for (int i = 0; i < static_cast<int>(ciphertext.size()); i++) {
    tmp = modular_exponentiation(ciphertext[i], d, n);
    result.push_back(tmp);
  }
  return result;
}

int main(int argc, char const *argv[]) {
  string tmp;
  string text;
  int size;
  cout << "Enter the size (in bits) of the key please : ";
  getline(cin, tmp);
  size = stoi(tmp);
  vector<Bignum> keys = key_generator(size);
  cout << "Enter text you want to encrypt please: ";
  getline(cin, text);
  cout << "Encrypting & decrypting the message ..." << endl;
  vector<Bignum> msg = encode(text, size);
  vector<Bignum> ciphertext = encrypt(msg, keys[0], keys[2]);
  vector<Bignum> plaintext = decrypt(ciphertext, keys[1], keys[2]);
  string decoded = decode(plaintext, size);
  cout << "Plaintext : " << decoded << endl;
  return 0;
}