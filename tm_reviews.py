from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition import NMF, LatentDirichletAllocation
from sklearn.datasets import fetch_20newsgroups
import numpy

n_samples = 1000
n_features = 2500
n_topics = 10
fname = "reviews.20k.txt"
#dataset = fetch_20newsgroups(shuffle=True, random_state=1,
#                             remove=('headers', 'footers', 'quotes'))
#data_samples = dataset.data[:n_samples]
with open(fname) as f:
    content = f.readlines()
dataset = [x.strip() for x in content] 
data_samples = dataset[:n_samples]
data_samples = [samp.split("\t")[-1] for samp in data_samples]

# Use tf-idf features for NMF.
print("Extracting tf-idf features for NMF...")
tfidf_vectorizer = TfidfVectorizer(max_df=0.95, min_df=2,
                                   max_features=n_features,
                                   stop_words='english')

tfidf = tfidf_vectorizer.fit_transform(data_samples)

# Fit the NMF model
print("Fitting the NMF model with tf-idf features, "
      "n_samples=%d and n_features=%d..."
      % (n_samples, n_features))

nmf = NMF(n_components=n_topics, random_state=1,
          alpha=.1, l1_ratio=.5).fit(tfidf)

tfidf_feature_names = tfidf_vectorizer.get_feature_names()

print("Writing tf-idf features and NMF model to file...")
vocab_file = open('vocabulary.txt', 'w')
p = 0
for item in tfidf_feature_names:
  try:
  	vocab_file.write("%s\n" % item)
  except:
  	p +=1
  	vocab_file.write("%s\n" % ("weird_symbol"+str(p)))

W = nmf.fit_transform(tfidf)
H = nmf.components_
numpy.savetxt("mat.csv", tfidf.toarray(), delimiter=",")
numpy.savetxt("W.csv", W, delimiter=",")
numpy.savetxt("H.csv", H, delimiter=",")

