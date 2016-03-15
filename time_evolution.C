void time_evolution () {

  int events[6]= {0,5000,10000,15000,20000,25000};
  int time[6]= {0,10,30,60,102,169};
  TGraph * g = new TGraph(6,events,time);
  g->Draw();
}
