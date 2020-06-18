from mcflow.changeProbs import main
import os

path = os.path.abspath(os.path.dirname(__file__))


def test():
    main(path=path, feeds=[os.path.join('change-probs', 'FAU-0')], verbosity=10, indep=(1, 2), debug=False)


if __name__ == '__main__':
    test()
