#include "utils.h"
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <version.h>
#include <cfloat>

static_assert(
    __cplusplus >= 201103L,
    "This code needs to be compiled with a c++11 compatible compiler");

class site {
public:
  string chr;
  int pos;
  string ref, alt;
  int idx;
  char keepNum;
  int type;
  vector<float> frq;

  explicit site(string &c, int p, string &r, string &a, int i)
      : site(c, p, r, a, i, 1) {}

  explicit site(string &c, int p, string &r, string &a, int i, char f)
      : chr(c), pos(p), ref(r), alt(a), idx(i), keepNum(f), type(1) {
    if ((ref == "A" || ref == "T" || ref == "G" || ref == "C") &&
        (alt == "A" || alt == "T" || alt == "G" || alt == "C"))
      type = 0;
  }

  bool strand(const string &r, const string &a) {
    if (ref == r && alt == a)
      return true;
    else
      return false;
  }
};

class site_map {
public:
  vector<site *> VS;

  // let's store all sites by chrom and then by position
  unordered_map<string, multimap<int, site *> > MS;

  site_map() {}

  void push(site *s) {
    VS.push_back(s);
    auto got = MS.find(s->chr);
    if (got == MS.end())
      got = MS.insert(make_pair(s->chr, multimap<int, site *>())).first;

    got->second.insert(pair<int, site *>(s->pos, s));
  }

  site *get(const string &chr, int pos, const string &ref, const string &alt) {

    // grab the correct map based on chromosome input
    auto chrMMapItr = MS.find(chr);
    if (chrMMapItr == MS.end())
      throw runtime_error("Could not find chromosome: " + chr);

    // grab the positions on the chromosome that match
    auto ret = chrMMapItr->second.equal_range(pos);

    for (auto it = ret.first; it != ret.second; ++it) {
      if (it->second->strand(ref, alt))
        return it->second;
    }
    return nullptr;
  }

  int size() { return VS.size(); }
};

void toLowerCase(string &str) {
  for (int c = 0; c < str.size(); c++)
    str[c] = tolower(str[c]);
}

int main(int argc, char **argv) {

  cout << "measureAggregateRsquared v" << VERSION_MAJOR << "." << VERSION_MINOR
       << endl;
  string buffer;
  vector<string> tok;

  char *validation = nullptr;
  char *imputed = nullptr;
  char *sample = nullptr;
  vector<char *> exclude;
  vector<char *> include;
  char *freq = nullptr;
  char *fmap = nullptr;
  char *output = nullptr;
  char *bin = nullptr;
  char keepNum = 0;
  bool noMapDump = false;
  bool discardFrqDoubles = false;
  bool discardMonomorphic = false;
  bool ignoreSecondFrqDouble = false;

  // reading arguments
  for (int a = 1; a < argc; a++) {
    if (strcmp(argv[a], "--validation") == 0) {
      validation = argv[a + 1];
      keepNum |= 1;
    };
    if (strcmp(argv[a], "--imputed") == 0) {
      imputed = argv[a + 1];
      keepNum |= 2;
    };
    if (strcmp(argv[a], "--sample") == 0) {
      sample = argv[a + 1];
    };
    if (strcmp(argv[a], "--exclude") == 0) {
      exclude.push_back(argv[a + 1]);
    };
    if (strcmp(argv[a], "--include") == 0) {
      include.push_back(argv[a + 1]);
    };
    if (strcmp(argv[a], "--freq") == 0) {
      freq = argv[a + 1];
      keepNum |= 4;
    };
    if (strcmp(argv[a], "--fmap") == 0) {
      fmap = argv[a + 1];
      keepNum |= 4;
    };
    if (strcmp(argv[a], "--output") == 0) {
      output = argv[a + 1];
    };
    if (strcmp(argv[a], "--bin") == 0) {
      bin = argv[a + 1];
    };
    if (strcmp(argv[a], "--no-map-dump") == 0) {
      noMapDump = true;
    }
    if (strcmp(argv[a], "--discard-fmap-doubles") == 0) {
      discardFrqDoubles = true;
    }
    if (strcmp(argv[a], "--ignore-second-fmap-double") == 0) {
      ignoreSecondFrqDouble = true;
    }
    if (strcmp(argv[a], "--discard-monomorphic") == 0) {
      discardMonomorphic = true;
    };
  }
  assert(keepNum == 7);

  // checking arguments
  cout << endl;
  if (validation != nullptr)
    cout << "Validation:\t" << validation << endl;
  else {
    cerr << "Argument --validation missing!" << endl;
    exit(1);
  }
  if (imputed != nullptr)
    cout << "Imputed:\t" << imputed << endl;
  else {
    cerr << "Argument --imputed missing!" << endl;
    exit(1);
  }
  if (sample != nullptr)
    cout << "Sample:\t\t" << sample << endl;
  else {
    cerr << "Argument --sample missing!" << endl;
    exit(1);
  }
  for (int e = 0; e < exclude.size(); e++)
    cout << "ExcludeSites:\t" << exclude[e] << endl;
  for (int e = 0; e < include.size(); e++)
    cout << "IncludeSites:\t" << include[e] << endl;

  if (freq != nullptr && fmap != nullptr) {
    cerr << "Please only specify --freq or --fmap";
    exit(1);
  } else if (freq != nullptr)
    cout << "Frequencies:\t" << freq << endl;
  else if (fmap != nullptr)
    cout << "Frequencies:\t" << fmap << endl;
  else {
    cerr << "Argument --freq or --fmap missing!" << endl;
    exit(1);
  }

  if (bin != nullptr)
    cout << "Bins:\t\t" << bin << endl;
  else {
    cerr << "Argument --bin missing!" << endl;
    exit(1);
  }
  if (output != nullptr)
    cout << "Output:\t\t" << output << endl;
  else {
    cerr << "Argument --output missing!" << endl;
    exit(1);
  }
  cout << endl;

  // Reading Exclusion
  unordered_set<pair<string, int>, boost::hash<std::pair<string, int> > > EXC;
  for (int e = 0; e < exclude.size(); e++) {
    cout << "Reading sites to exclude in [" << exclude[e] << "]" << endl;
    ifile fde(exclude[e]);
    while (getline(fde, buffer, '\n')) {
      tok = sutils::tokenize(buffer, "\t");
      if (tok.size() != 2)
        throw runtime_error(
            "--exclude file " + string(exclude[e]) +
            " does not contain sites in 'chrom\tposition' format");
      EXC.insert(pair<string, int>(tok[0], stoi(tok[1])));
    }
    fde.close();
    cout << "  * #sites=" << EXC.size() << endl;
  }

  // Reading Inclusion
  unordered_set<pair<string, int>, boost::hash<std::pair<string, int> > > INC;
  for (int e = 0; e < include.size(); e++) {
    cout << "Reading sites to include in [" << include[e] << "]" << endl;
    ifile fdi(include[e]);
    while (getline(fdi, buffer, '\n')) {
      tok = sutils::tokenize(buffer, "\t");
      if (tok.size() != 2)
        throw runtime_error(
            "--include file " + string(include[e]) +
            " does not contain sites in 'chrom\tposition' format");
      INC.insert(pair<string, int>(tok[0], stoi(tok[1])));
    }
    fdi.close();
    cout << "  * #sites=" << INC.size() << endl;
  }

  // Reading Samples
  int idx_pop = -1;
  map<string, vector<int> > S;
  cout << "Reading samples description in [" << sample << "]" << endl;
  ifile fds(sample);
  getline(fds, buffer, '\n');
  tok = sutils::tokenize(buffer, " ");
  for (int t = 0; t < tok.size() && idx_pop < 0; t++) {
    toLowerCase(tok[t]);
    if (tok[t] == "pop" || tok[t] == "population")
      idx_pop = t;
  }
  if (idx_pop < 0) {
    cerr << "You must have a pop column in the sample file" << endl;
    exit(1);
  } else
    cout << "  * population column: " << idx_pop << endl;
  getline(fds, buffer, '\n');
  int n_ind = 0;
  while (getline(fds, buffer, '\n')) {
    tok = sutils::tokenize(buffer, " ");
    map<string, vector<int> >::iterator itS = S.find(tok[idx_pop]);
    if (itS != S.end())
      itS->second.push_back(n_ind);
    else
      S.insert(pair<string, vector<int> >(tok[idx_pop], vector<int>(1, n_ind)));
    n_ind++;
  }
  cout << "  * #pop=" << S.size() << endl;
  cout << "  * #ind=" << n_ind << endl;
  cout << "  * effectif=";
  for (map<string, vector<int> >::iterator itS = S.begin(); itS != S.end();
       itS++)
    cout << " [" << itS->first << ", " << itS->second.size() << "]";
  cout << endl;

  // Reading Validation data
  int i_site = 0;
  site_map SM;
  vector<vector<float> > DV;
  cout << "Reading validation data in [" << validation << "]" << endl;
  ifile fdv(validation);
  while (getline(fdv, buffer, '\n')) {
    tok = sutils::tokenize(buffer, " ");
    assert(tok.size() == (n_ind * 3 + 5));
    int pos = atoi(tok[2].c_str());
    string chr = tok[0];

    auto itE = EXC.find(make_pair(chr, pos));
    auto itI = INC.find(make_pair(chr, pos));

    // don't bother storing sites that we don't want to keep
    char flag = 1;
    if ((!INC.empty() && itI == INC.end()) ||
        (!EXC.empty() && itE != EXC.end())) {
      if (fmap != nullptr)
        continue;
      else
        flag = 0;
    }

    site *s = new site(chr, pos, tok[3], tok[4], SM.size(), flag);
    SM.push(s);
    DV.push_back(vector<float>(n_ind, -1.0));
    for (int i = 5; i < tok.size(); i += 3) {
      double g0 = atof(tok[i + 0].c_str());
      double g1 = atof(tok[i + 1].c_str());
      double g2 = atof(tok[i + 2].c_str());
      if (g0 >= 0.9 || g1 >= 0.9 || g2 >= 0.9)
        DV.back()[(i - 5) / 3] = g1 + 2 * g2;
    }
    if (i_site % 1000 == 0) {
      cout << "\r" << i_site;
    }
    i_site++;
  }
  cout << "\r" << i_site << endl;
  cout << "  * #site=" << SM.size() << endl;
  fdv.close();

  // Reading Frequencies
  vector<string> P;
  P.reserve(SM.size());
  if (freq != nullptr) {
    cout << "Reading frequencies in [" << freq << "]" << endl;

    ifile fdf(freq);
    getline(fdf, buffer, '\n');
    tok = sutils::tokenize(buffer, " ");
    for (int t = 0; t < tok.size(); t++)
      P.push_back(tok[t]);
    cout << "  * #pop to be analysed: " << P.size() << endl;
    for (int s = 0; s < SM.size(); s++) {
      if (!getline(fdf, buffer, '\n')) {
        cerr << "Frequency file must contain same number of lines than "
                "validation data + 1 (not enough lines)" << endl;
        exit(1);
      }
      tok = sutils::tokenize(buffer, " ");
      for (int t = 0; t < tok.size(); t++)
        SM.VS[s]->frq.push_back(atof(tok[t].c_str()));
      SM.VS[s]->keepNum |= 4;
    }
    if (getline(fdf, buffer, '\n')) {
      cerr
          << "Frequency file must contain same number of lines than validation "
             "data + 1 (too much lines)" << endl;
      exit(1);
    }
    fdf.close();
  } else if (fmap != nullptr) {
    cout << "Reading frequencies in [" << fmap << "]" << endl;

    ifile fdf(fmap);
    getline(fdf, buffer, '\n');
    tok = sutils::tokenize(buffer, "\t");
    for (int t = 4; t < tok.size(); t++)
      P.push_back(tok[t]);
    cout << "  * #pop to be analysed: " << P.size() << endl;
    unsigned nFreqFound = 0;
    unsigned nFreq = 0;
    while (getline(fdf, buffer, '\n')) {
      tok = sutils::tokenize(buffer, "\t");
      auto s = SM.get(tok[0], stoul(tok[1]), tok[2], tok[3]);
      if (s) {

        // if the freq file contains sites more than once, then throw them out
        // or die
        if (!s->frq.empty()) {
          if (discardFrqDoubles) {
              s->keepNum = 0;
            continue;
          } else if (ignoreSecondFrqDouble)
            continue;
          else
            throw runtime_error(
                "fmap file contains duplicate freq lines at site: " + buffer);
        }
        for (int t = 4; t < tok.size(); t++)
          s->frq.push_back(stof(tok[t]));
        ++nFreqFound;
        s->keepNum |= 4;
      }
      ++nFreq;
    }
    cout << "  * #frequency sites found: " << nFreq << endl;
    cout << "  * #frequency sites kept: " << nFreqFound << endl;
    fdf.close();
  } else
    assert(false);

  // Reading Imputation Data
  cout << "Reading imputation data in [" << imputed << "]" << endl;
  i_site = 0;
  int n_excluded = 0;
  int n_kept = 0;
  vector<vector<float> > DI = DV;
  ifile fdi(imputed);
  while (getline(fdi, buffer, '\n')) {
    tok = sutils::tokenize(buffer, " ");
    assert(tok.size() == (n_ind * 3 + 5));

    site *s = SM.get(tok[0], atoi(tok[2].c_str()), tok[3], tok[4]);
    if (s != nullptr) {
      for (int i = 5; i < tok.size(); i += 3)
        DI[s->idx][(i - 5) / 3] =
            atof(tok[i + 1].c_str()) + 2 * atof(tok[i + 2].c_str());
      s->keepNum |= 2;
      ++n_kept;
    } else
      ++n_excluded;

    if (i_site % 1000 == 0)
      cout << "\r" << i_site;

    i_site++;
  }

  cout << "\r" << i_site << endl;
  cout << "  * #site=" << i_site << endl;
  cout << "  * #included=" << n_kept << endl;
  cout << "  * #excluded=" << n_excluded << endl;
  fdi.close();

  // Reading Bins
  vector<double> B;
  cout << "Reading frequency bins in [" << bin << "]" << endl;
  ifile fdb(bin);
  while (getline(fdb, buffer, '\n'))
    B.push_back(atof(buffer.c_str()));
  fdb.close();
  int n_bins = B.size();
  if (n_bins < 2)
    throw runtime_error(
        "Please specify at least two boundaries in the bins file");

  assert(B[0] <= 1 && B[0] >= 0);
  if (!discardMonomorphic && B[0] == 0)
    B[0] -= DBL_EPSILON;

  cout << "  * #bins=" << n_bins - 1 << endl;

  cout << "Summarizing data" << endl;
  int n_snp = 0;
  int n_imp = 0;
  int n_cmp = 0;
  for (int s = 0; s < SM.size(); s++) {
    if (SM.VS[s]->keepNum == keepNum) {
      if (SM.VS[s]->type == 0)
        n_snp++;
      else
        n_cmp++;
      n_imp++;
    }
  }
  cout << "  * #validation sites=" << SM.size() << endl;
  cout << "  * #at which to compare r-square=" << n_imp << endl;
  cout << "  * #snp=" << n_snp << endl;
  cout << "  * #complex=" << n_cmp << endl;

  if (n_snp + n_cmp == 0)
    throw runtime_error(
        "No snps or complex variants were found to do an analysis on");

  set<int> countedSites;

  // Measuring Aggregate Rsquared per population x type x bin
  const int numTypes = 3;
  vector<vector<vector<double> > > A = vector<vector<vector<double> > >(
      P.size(),
      vector<vector<double> >(numTypes, vector<double>(n_bins - 1, -1.0)));
  vector<vector<vector<double> > > D = vector<vector<vector<double> > >(
      P.size(),
      vector<vector<double> >(numTypes, vector<double>(n_bins - 1, -1.0)));
  vector<vector<vector<double> > > F = vector<vector<vector<double> > >(
      P.size(),
      vector<vector<double> >(numTypes, vector<double>(n_bins - 1, -1.0)));

  // aggregate allele frequency of validation data in each bin
  vector<vector<vector<double> > > freqs_validation = vector<vector<vector<double> > >(
      P.size(),
      vector<vector<double> >(numTypes, vector<double>(n_bins - 1, -1.0)));

  // aggregate allele frequency of imputed data in each bin
  vector<vector<vector<double> > > freqs_imputed = vector<vector<vector<double> > >(
      P.size(),
      vector<vector<double> >(numTypes, vector<double>(n_bins - 1, -1.0)));
  
  
  for (int p = 0; p < P.size(); p++) {
    map<string, vector<int> >::iterator itS = S.find(P[p]);
    if (itS != S.end()) {
      vector<int> I = itS->second;
      for (int t = 0; t < numTypes; t++) { // t=2 means both snp and complex
        if ((t == 0 && n_snp > 0) || (t == 1 && n_cmp > 0) ||
            (t == 2 && n_snp + n_cmp > 0)) {
          for (int b = 1; b < n_bins; b++) {
            cout << "Rsquared\t[pop=" << P[p] << "]";
            cout << "\t[type=";
            switch (t) {
            case 0:
              cout << "SNPs]";
              break;
            case 1:
              cout << "COMPLEXs]";
              break;
            case 2:
              cout << "ALL]";
              break;
            }
            cout << "\t[bin=" << sutils::double2str(B[b], 3) << "]";

            // Calculate Means
            double mean_frq = 0.0;
            int count_frq = 0;
            double mean_ref = 0.0;
            double mean_exp = 0.0;
            int mean_cnt = 0;

            for (int s = 0; s < SM.size(); s++) {
              if (SM.VS[s]->keepNum == keepNum &&
                  (t == 2 || SM.VS[s]->type == t) &&
                  SM.VS[s]->frq[p] > B[b - 1] && SM.VS[s]->frq[p] <= B[b]) {
                countedSites.insert(s);
                mean_frq += SM.VS[s]->frq[p];
                count_frq++;
                for (int i = 0; i < I.size(); i++) {
                  if (DV[s][I[i]] >= 0) {
                    mean_ref += DV[s][I[i]];
                    mean_exp += DI[s][I[i]];
                    mean_cnt++;
                  }
                }
              }
            }

            if (mean_cnt == 0) {
              cerr << "Calculation aborted, number of validation genotypes = 0"
                   << endl;
            } else {
              cout << "\t[geno=" << mean_cnt << "]";
              mean_ref /= mean_cnt;
              mean_exp /= mean_cnt;
              cout << "\t[me_r=" << sutils::double2str(mean_ref, 4) << "]";
              cout << "\t[me_e=" << sutils::double2str(mean_exp, 4) << "]";

              // Calculate Standard Deviations
              double std_ref = 0.0;
              double std_exp = 0.0;

              for (int s = 0; s < SM.size(); s++) {
                if (SM.VS[s]->keepNum == keepNum &&
                    (t == 2 || SM.VS[s]->type == t) &&
                    SM.VS[s]->frq[p] > B[b - 1] && SM.VS[s]->frq[p] <= B[b]) {
                  for (int i = 0; i < I.size(); i++) {
                    if (DV[s][I[i]] >= 0) {
                      std_ref +=
                          (DV[s][I[i]] - mean_ref) * (DV[s][I[i]] - mean_ref);
                      std_exp +=
                          (DI[s][I[i]] - mean_exp) * (DI[s][I[i]] - mean_exp);
                    }
                  }
                }
              }
              std_ref = sqrt(std_ref / (mean_cnt - 1));
              std_exp = sqrt(std_exp / (mean_cnt - 1));
              cout << "\t[sd_r=" << sutils::double2str(std_ref, 4) << "]";
              cout << "\t[sd_e=" << sutils::double2str(std_exp, 4) << "]";

              // Calculate Rsquared
              double sum = 0.0;
              for (int s = 0; s < SM.size(); s++) {
                if (SM.VS[s]->keepNum == keepNum &&
                    (t == 2 || SM.VS[s]->type == t) &&
                    SM.VS[s]->frq[p] > B[b - 1] && SM.VS[s]->frq[p] <= B[b]) {
                  for (int i = 0; i < I.size(); i++) {
                    if (DV[s][I[i]] >= 0) {
                      sum += ((DV[s][I[i]] - mean_ref) / std_ref) *
                             ((DI[s][I[i]] - mean_exp) / std_exp);
                    }
                  }
                }
              }
              sum /= (mean_cnt - 1);
              A[p][t][b - 1] = sum * sum;
              D[p][t][b - 1] = mean_frq / count_frq;
              F[p][t][b - 1] = count_frq;
              freqs_validation[p][t][b - 1] =
                  mean_ref /
                  2; // average dose over 2 == validation allele frequency
              freqs_imputed[p][t][b - 1] =
                  mean_exp /
                  2; // average dose over 2 == validation allele frequency

              cout << "\t[site=" << count_frq << "]";
              cout << "\t[frq=" << sutils::double2str(mean_frq / count_frq, 4)
                   << "]";
              cout << "\t[rsq=" << sutils::double2str(A[p][t][b - 1], 4) << "]"
                   << endl;
              ;
            }
          }
        }
      }
    }
  }

  // Writing Results
  for (int p = 0; p < P.size(); p++) {
    map<string, vector<int> >::iterator itS = S.find(P[p]);
    if (itS != S.end()) {
      vector<int> I = itS->second;
      for (int t = 0; t < numTypes; t++) {
        if ((t == 0 && n_snp > 0) || (t == 1 && n_cmp > 0) ||
            (t == 2 && n_snp + n_cmp > 0)) {
          string filename = string(output) + "." + P[p];
          switch (t) {
          case 0:
            filename += ".snps";
            break;
          case 1:
            filename += ".complexs";
            break;
          case 2:
            filename += ".all";
            break;
          }

          cout << "Writing results in [" << filename << "]" << endl;
          ofile fdo(filename.c_str());
          fdo << "Bin_frequency\tr_square\tnum_genotypes\tfreq_validation\tfreq_imputation" << endl;
          for (int b = 1; b < n_bins; b++)
            fdo << D[p][t][b - 1] << "\t" << A[p][t][b - 1] << "\t"
                << F[p][t][b - 1] << "\t" << freqs_validation[p][t][b - 1] << "\t"
                << freqs_imputed[p][t][b - 1] << endl;
          fdo.close();

          // write map of sites used
          if (!noMapDump) {
            string mapFile = filename + ".map";
            cout << "Writing map of all compared sites to [" << mapFile << "]"
                 << endl;
            ofile fdmap(mapFile.c_str());
            for (auto sNum : countedSites) {
              auto s = SM.VS[sNum];
              assert(s);
              fdmap << s->chr << "\t" << to_string(s->pos) << "\t" << s->ref
                    << "\t" << s->alt << "\n";
            }
            fdmap.close();
          }
        }
      }
    }
  }
}
