
#ifndef POWER_H
#define POWER_H

class Power{
private:
  int K;

public:
  Power(){K = 1;}; // default constructor
  Power(int k);
  double xToTheK(double x);
};

#endif /* POWER_H */
