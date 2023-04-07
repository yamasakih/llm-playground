import json
import openai
from typing import Dict, Any, List


def calculate_embedding_mean(json_data: Dict[Any, Any]) -> List[float]:
    SMILES = json_data["SMILES"]
    model = "text-embedding-ada-002"

    embeddings = []
    for smiles in SMILES:  # マニュアルで smile -> smiles に変更
        response = openai.Embedding.create(input=smiles, model=model)
        embedding = response["data"][0]["embedding"]
        embeddings.append(embedding)

    # 要素ごとの平均を計算
    embedding_mean = [sum(x) / len(x) for x in zip(*embeddings)]

    return embedding_mean


# ファイルを読み込む
with open("resources/CBLB_inhibitors_vsF.json", "r") as f:
    data = json.load(f)

# 各要素に対して calculate_embedding_mean を実行し、embedding_mean を追加する
for item in data:
    item["embedding_mean"] = calculate_embedding_mean(item)

# 結果をファイルに保存する
with open("resources/CBLB_inhibitors_vsF_added_embedding_mean.json", "w") as f:
    json.dump(data, f)
