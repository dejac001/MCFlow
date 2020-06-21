from changeProbs import main
import os

path = os.path.abspath(os.path.dirname(__file__))


def test():
    main(path=os.path.join(path, 'change-probs'), feeds=[os.path.join('FAU-0')], verbosity=10, indep=(1, 2),
         debug=False)
    main(path=os.path.join(path, 'change-probs'), feeds=[os.path.join('FAU-0')], verbosity=10, indep=(1, 2),
         debug=False)


if __name__ == '__main__':
    test()
