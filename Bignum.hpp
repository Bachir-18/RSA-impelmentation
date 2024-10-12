#ifndef BIGNUM_HPP
#define BIGNUM_HPP

#include <string>
#include <vector>

class Bignum {
public:
  Bignum(int = 0);
  Bignum(unsigned);
  Bignum(std::string s);
  Bignum(Bignum const &) = default;
  Bignum(Bignum &&) = default;

  ~Bignum() = default;

  Bignum &operator=(Bignum const &) = default;
  Bignum &operator=(Bignum &&) = default;

  Bignum &operator+=(Bignum const &x) {
    *this = *this + x;
    return *this;
  };
  Bignum &operator-=(Bignum const &x) {
    *this = *this - x;
    return *this;
  };
  Bignum &operator*=(Bignum const &x) {
    *this = *this * x;
    return *this;
  };
  Bignum &operator/=(Bignum const &x) {
    *this = *this / x;
    return *this;
  };
  Bignum &operator<<=(unsigned n);
  Bignum &operator>>=(unsigned n);

  uint32_t &operator[](unsigned i) { return tab[i]; };
  uint32_t operator[](unsigned i) const { return tab[i]; };

  size_t size() const { return tab.size(); };

  void printHex() const;
  void resize();

private:
  std::vector<uint32_t> tab;
  bool isPositive;

  friend bool operator<(Bignum const &, Bignum const &);
  friend bool operator==(Bignum const &, Bignum const &);
  friend bool operator<=(Bignum const &, Bignum const &);
  friend bool operator>(Bignum const &x, Bignum const &y) { return !(x <= y); };
  friend bool operator>=(Bignum const &x, Bignum const &y) { return !(x < y); };
  friend bool operator!=(Bignum const &x, Bignum const &y) {
    return !(x == y);
  };

  friend std::ostream &operator<<(std::ostream &, Bignum const &);

  friend Bignum operator<<(Bignum const &x, unsigned n) {
    auto y = x;
    return y <<= n;
  };
  friend Bignum operator>>(Bignum const &x, unsigned n) {
    auto y = x;
    return y >>= n;
  };

  friend Bignum operator+(Bignum const &, Bignum const &);
  friend Bignum operator-(Bignum const &, Bignum const &);
  friend Bignum operator*(Bignum const &, Bignum const &);
  friend Bignum operator/(Bignum const &x, Bignum const &y);
  friend Bignum operator%(Bignum const &, Bignum const &);

  friend std::pair<Bignum, Bignum> division(Bignum const &, Bignum const &);

  friend Bignum inverseMod(Bignum const &, Bignum const &);

  friend bool compareAbs(Bignum const &x, Bignum const &y);
  friend Bignum SubtractX_Y(Bignum const &x, Bignum const &y);
  friend Bignum addSameSign(Bignum const &x, Bignum const &y);
  friend Bignum prime_number_generator(unsigned int n);
  friend Bignum odd_bignum_generator(unsigned int n);
  friend Bignum modular_exponentiation(const Bignum &a, const Bignum &b,
                                       const Bignum &n);
  friend Bignum generateRandomPrime(int bits);
  friend std::string decode(std::vector<Bignum> encodedText, int n);
};

#endif
