import pandas as pd  
import numpy as np  
import matplotlib.pyplot as plt  
import seaborn as sb 
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LinearRegression
from sklearn import metrics

df = pd.read_csv("all_log_results.csv", index_col='gid')
df.shape
print(df.describe())

x_cols = ["avg_exc", "avg_inh", "max_exc", "max_inh", "num_exc", "num_inh"]

X = df[x_cols].values
y = df['FR'].values

# plt.figure()
# sb.distplot(df['FR'])
# plt.show()

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

regressor = LinearRegression()
regressor.fit(X_train, y_train)

coeff_df = pd.DataFrame(regressor.coef_, x_cols, columns=['Coefficient'])
print(coeff_df)

y_pred = regressor.predict(X_test)

comp = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
print(comp.head(25))

plt.figure()
plt.hist(y_pred, label="preds", alpha = 0.5)
plt.hist(y_test, label='trues', alpha = 0.5)
plt.legend()
plt.show()

plt.figure()
comp.head(25).plot(kind='bar',figsize=(10,8))
plt.grid(which='major', linestyle='-', linewidth='0.5', color='green')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.show()

print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))  
print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))  
print('Root Mean Squared Error:', np.sqrt(metrics.mean_squared_error(y_test, y_pred)))
print("R2 Score:", metrics.r2_score(y_test, y_pred))
