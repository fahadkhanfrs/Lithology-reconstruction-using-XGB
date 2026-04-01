from sklearn.model_selection import train_test_split

def preprocess(df):
    cols = ["Depth", "GR", "SP", "RILD", "CNLS", "DPOR", "RHOB"]
    df = df[cols].dropna()

    # Sort by depth (very important)
    df = df.sort_values(by="Depth")

    # Split by depth (not random)
    split = int(0.8 * len(df))

    train_df = df.iloc[:split]
    test_df = df.iloc[split:]

    X_train = train_df.drop(["RHOB", "Depth"], axis=1)
    y_train = train_df["RHOB"]

    X_test = test_df.drop(["RHOB", "Depth"], axis=1)
    y_test = test_df["RHOB"]

    return X_train, X_test, y_train, y_test