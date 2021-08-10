import math
import copy
from functools import wraps
import time


def stop_watch(func):
    @wraps(func)
    def wrapper(*args, **kargs):
        start = time.time()
        result = func(*args, **kargs)
        process_time = time.time() - start
        print(f"{func.__name__}の実行時間:{process_time}[sec]")
        return result
    return wrapper


def score(a, b):
    """入力された2塩基間のマッチ/アンマッチスコアを返す

    Args:
        a (str):　スコアを計算する塩基(長さ1の文字)
        b (str):　スコアを計算する塩基(長さ1の文字)

    Returns:
        int: 塩基a, b間のマッチ/アンマッチスコア
    """
    return 2 if a == b else -1


@stop_watch
def alignment(seq_1, seq_2, gap_start, gap_extention):
    """DP法を用いて`seq_1`と`seq_2`の配列間でアラインメントを行う

    Args:
        seq_1 (str): アラインメントしたいDNA配列
        seq_2 (str): アラインメントしたいDNA配列
        gap_start (int): Affine Gapにおけるギャップ開始ペナルティ
        gap_extention (int): Affine Gapにおけるギャップ伸長ペナルティ
    """

    X = len(seq_1)  # 配列1の長さ
    Y = len(seq_2)  # 配列2の長さ

    # 累積スコアを記録する二次元配列を用意する。負の無限大で埋めておく
    dp = [[-math.inf]*(Y+1) for _ in range(X+1)]

    # 累積スコアがどの方向から計算されたものかを記録する二次元配列
    directions = [[0]*(Y+1) for _ in range(X+1)]
    # 0 で左から
    # 1 で上から
    # 2 で左上から計算された事を表す

    # 各配列の初期化
    for i in range(X+1):
        dp[i][0] = 0          # 累積スコア配列の上辺を 0 で初期化する(out gap penalty = 0)
        directions[i][0] = 0  # 上辺では左方向(0)からのみ移動できるので 0 で初期化する

    for j in range(Y+1):
        dp[0][j] = 0          # 累積スコア配列の左辺を 0 で初期化する(out gap penalty = 0)
        directions[0][j] = 1  # 左辺では上方向(1)からのみ移動できるので 1 で初期化する

    # ギャップが連続しているかを記憶する配列(affine gapの計算用)を用意する
    # 配列のサイズと初期条件が`dp`と同一であるためcopyを用いた
    gap_auxiliary_X = copy.deepcopy(dp)  # X方向(横方向)の補助配列
    gap_auxiliary_Y = copy.deepcopy(dp)  # Y方向(縦方向)の補助配列

    # Matrix Fill step
    for i in range(1, X+1):
        for j in range(1, Y+1):

            # ギャップの連続判定用補助配列を漸化式に従って計算する
            gap_auxiliary_X[i][j] = max(dp[i-1][j]+gap_start,
                                        gap_auxiliary_X[i-1][j]+gap_extention)
            gap_auxiliary_Y[i][j] = max(dp[i][j-1]+gap_start,
                                        gap_auxiliary_Y[i][j-1]+gap_extention)

            # 累積スコアを漸化式に従って計算する
            # 漸化式に現れる最大値の計算のため、配列 scores を用意して max() を利用する
            if i == X:
                # 右辺では上方向からの移動が out gap を表すので、ペナルティを 0 とする
                scores = [gap_auxiliary_X[i][j],
                          dp[i][j-1],
                          dp[i-1][j-1]+score(seq_1[i-1], seq_2[j-1])]
            elif j == Y:
                # 下辺では左方向からの移動が out gap を表すので、ペナルティを 0 とする
                scores = [dp[i-1][j],
                          gap_auxiliary_Y[i][j],
                          dp[i-1][j-1]+score(seq_1[i-1], seq_2[j-1])]
            else:
                scores = [gap_auxiliary_X[i][j],
                          gap_auxiliary_Y[i][j],
                          dp[i-1][j-1]+score(seq_1[i-1], seq_2[j-1])]

            # [(0, 左から来た場合のスコア),
            #  (1, 上から来た場合のスコア),
            #  (2, 左上から来た場合のスコア)]
            # の形の配列からスコアが最大となるものを探し、
            # 方向を表す整数と累積スコアをそれぞれ配列に記録する
            direction, max_score = max(enumerate(scores), key=lambda x: x[1])
            dp[i][j] = max_score
            directions[i][j] = direction

    # 累積スコア配列の終点がGlobal Alignmentのスコアとなる
    global_score = dp[X][Y]
    print(f"global score:{global_score}")

    # Trace Back step

    pairs = []  # アラインメントの結果ペアと判定した文字列の配列を格納する配列

    # directions[X][Y]からdirections[0][0]まで累積スコア計算時の方向を逆にたどる
    i = X
    j = Y
    while i != 0 or j != 0:
        if directions[i][j] == 0:
            # 左から計算された場合は配列2側にギャップを挿入する
            pairs.append([seq_1[i-1], " ", "-"])
            i -= 1
        elif directions[i][j] == 1:
            # 上から計算された場合は配列1側にギャップを挿入する
            pairs.append(["-", " ", seq_2[j-1]])
            j -= 1
        else:
            # 左上から計算された場合はその位置にある両塩基をペアとする
            if seq_1[i-1] == seq_2[j-1]:
                pairs.append([seq_1[i-1], "|", seq_2[j-1]])
            else:
                pairs.append([seq_1[i-1], " ", seq_2[j-1]])
            i -= 1
            j -= 1

    pairs.reverse()  # 逆向きに辿っていたので、正しい順番になるよう配列を逆転する

    # pairs を見やすく整列し、文字列にする
    align_text = "\n".join(["".join(seq_list) for seq_list in [
                           list(seq_tuple) for seq_tuple in zip(*pairs)]])

    print(align_text)  # アラインメント結果を描画する
    return


# 大腸菌
seq_1 = "ATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGCACAAAGAGGGGGACCTTAGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCAACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCNGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGANNTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACTTCGGGAGGGCG"

# ペスト菌
seq_2 = "ATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGCAGCGGGAAGTAGTTTACTACTTTGCCGGCGAGCGGCGGACGGGTGAGTAATGTCTGGGGATCTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATGACCTCGCAAGAGCAAAGTGGGGGACCTTAGGGCCTCACGCCATCGGATGAACCCAGATGGGATTAGCTAGTAGGTGGGGTAATGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTTCAGCGAGGAGGAAGGGGTTGAGTTTAATACGCTCAATCATTGACGTTACTCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGCGCTTAACGTGGGAACTGCATTTGAAACTGGCAAGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTACTCTTGACATCCACAGAATTTGGCAGAGATGCTAAAGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCACGTAATGGTGGGAACTCAAGGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGAGTAGGGCTACACACGTGCTACAATGGCAGATACAAAGTGAAGCGAACTCGCGAGAGCCAGCGGACCACATAAAGTCTGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTAGATCAGAATGCTACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGG"

alignment(seq_1, seq_2, -3, -1)
