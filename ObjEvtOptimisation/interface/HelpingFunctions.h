class ptGreater {
  public:
  template <typename T> bool operator () (const T* i, const T* j) {
    return (i->pt() > j->pt());
  }
};
