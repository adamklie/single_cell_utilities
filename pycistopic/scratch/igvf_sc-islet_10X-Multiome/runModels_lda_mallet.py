import pickle
import sys
import argparse
import os
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *

def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Run topic models.",)
    parser.add_argument('--inputcisTopic_obj', '-i', type=str, required=True,
                        help='Path to cisTopic object pickle file.')
    parser.add_argument('--output', '-o', type=str, required=False, default='models.pkl',
                        help='Path to save final model list.')
    parser.add_argument('--n_topics', '-nt', type=str, required=True, nargs='+',
                        help='Txt file containing selected topic id.')
    parser.add_argument('--n_cpu', '-c', type=int, required=True,
                        help = 'Number of cores')
    parser.add_argument('--n_iter', '-it', type=int, required=False, default=150,
                        help = 'Number of iterations')
    parser.add_argument('--alpha', '-a', type=int, required=False,  default=50,
                        help='Alpha value')
    parser.add_argument('--alpha_by_topic', '-abt', type=str, required=False, default=True,
                        help = 'Whether the alpha value should by divided by the number of topics')
    parser.add_argument('--eta', '-e', type=float, required=False, default=0.1,
                        help='Eta value.')
    parser.add_argument('--eta_by_topic', '-ebt', type=str, required=False, default=False,
                        help = 'Whether the eta value should by divided by the number of topics')
    parser.add_argument('--save_path', '-sp', type=str, required=False,
                        default=None, help='Whether intermediate models should be saved')
    parser.add_argument('--seed', '-s', type=int, required=False,
                        default=555, help='Seed for ensuring reproducibility')
    parser.add_argument('--temp_dir', '-td', type=str, required=False,
                        default=None, help='Path to TMP directory')
    parser.add_argument('--reuse_corpus', '-rc', type=str, required=False, default=False,
                        help = 'Whether to reuse the corpus from Mallet')
    parser.add_argument('--mallet_binary', '-mb', type=str, required=False, default='/cellar/users/aklie/opt/Mallet-202108/bin/mallet', 
                        help = 'Path to Mallet binary')
    return parser

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    """
    The main executable function
    """
    
    parser = make_argument_parser()
    args = parser.parse_args()
    filename = args.inputcisTopic_obj
    infile = open(filename, 'rb')
    cisTopic_obj = pickle.load(infile)
    infile.close()
    print('Input cisTopic_object:', filename)
    print(f'{cisTopic_obj.project} cisTopic object loaded') 
    print(f'cistopic_obj_load.fragment_matrix.shape: {cisTopic_obj.fragment_matrix.shape}')

    output = args.output
    print('Output file:', output)
    
    n_topics = args.n_topics
    print('Number of topics:', n_topics)
    n_topics = list(map(lambda x: int(x.split(',')[0]), n_topics))
    print('Number of topics:', n_topics)
    
    alpha=args.alpha
    print('Alpha:', alpha)
    
    alpha_by_topic=str2bool(args.alpha_by_topic)
    print('Divide alpha by the number of topics:', alpha_by_topic)
    
    eta=args.eta
    print('Eta:', eta)
    
    eta_by_topic=str2bool(args.eta_by_topic)
    print('Divide eta by the number of topics:', eta_by_topic)
    
    n_iter=args.n_iter
    print('Number of iterations:', n_iter)
    
    n_cpu=args.n_cpu
    print('Number of cores:', n_cpu)
    
    save_path=args.save_path
    print('Path to save intermediate files:', save_path)
    if save_path == 'None':
        save_path = None
    
    random_state=args.seed
    print('Seed:', random_state)
    
    temp_dir=args.temp_dir
    if temp_dir[-1] != '/':
        temp_dir = temp_dir + '/'
    print('Path to TMP dir:', temp_dir)
    
    reuse_corpus=args.reuse_corpus
    print('Reuse Mallet corpus:', reuse_corpus)

    mallet_binary=args.mallet_binary
    print('Path to Mallet binary:', mallet_binary)
    
    # Run models
    print('Running models')
    print('--------------')

    os.environ['MALLET_MEMORY'] = '100G'
    
    models=run_cgs_models_mallet(path_to_mallet_binary=mallet_binary, 
                                 cistopic_obj=cisTopic_obj,
                                 n_topics=n_topics,
                                 n_cpu=n_cpu,
                                 n_iter=n_iter,
                                 random_state=random_state,
                                 alpha=alpha,
                                 alpha_by_topic=alpha_by_topic,
                                 eta=eta,
                                 eta_by_topic=eta_by_topic,
                                 save_path=save_path,
                                 top_topics_coh=5,
                                 tmp_path=temp_dir,
                                 reuse_corpus=reuse_corpus)

    # Save
    if len(n_topics) == 1 and save_path is not None:
        print(f"Not saving to {output} because only one model was run and save_path was specified")
    else:
        print(f"Saving to {output}")
        with open(output, 'wb') as f:
            pickle.dump(models, f)

if __name__ == "__main__":
    main()
    