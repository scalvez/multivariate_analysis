
int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double Poisson_probability(int n_obs, int n_exp)
{
  return pow(n_exp,n_obs)*exp(-1*n_exp)/factorial(n_obs);
}
